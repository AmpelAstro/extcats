from pathlib import Path

import pandas as pd
from healpy import ang2pix
import pymongo
import pytest
from extcats import testdata
from extcats.CatalogPusher import CatalogPusher
import subprocess
import json


@pytest.fixture(scope="session")
def data_dir():
    return Path(__file__).parent / "test-data"


@pytest.fixture(scope="session")
def mongo_db(pytestconfig):
    """
    Bring up an instance of mongo in Docker
    """
    # if not pytestconfig.getoption("--integration"):
    #     raise pytest.skip("integration tests require --integration flag")
    try:
        container_id = (
            subprocess.check_output(
                ["docker", "run", "-d", "--rm", "-P", "mongo:4.0"],
            )
            .strip()
            .decode()
        )
    except FileNotFoundError:
        raise pytest.skip("integration test requires docker")
    try:
        port = json.loads(
            subprocess.check_output(
                ["docker", "inspect", container_id],
            ).decode()
        )[0]["NetworkSettings"]["Ports"]["27017/tcp"][0]["HostPort"]
        yield f"mongodb://localhost:{port}"
    finally:
        subprocess.check_output(
            ["docker", "stop", container_id],
        )


@pytest.fixture(scope="session")
def mongo_client(mongo_db):
    return pymongo.MongoClient(mongo_db)


@pytest.fixture(scope="session")
def milliquas(data_dir, mongo_client):
    # obtained with:
    # wget https://heasarc.gsfc.nasa.gov/FTP/heasarc/dbase/tdat_files/heasarc_milliquas.tdat.gz
    # gzcat heasarc_milliquas.tdat.gz| head -n 100 > heasarc_milliquas_sample.tdat
    # echo '<END>' >> heasarc_milliquas_sample.tdat
    mqp = CatalogPusher(
        catalog_name="milliquas",
        data_source=str(data_dir / "milliquas"),
        file_type="tdat",
    )
    columns = [
        "name",
        "ra",
        "dec",
        "lii",
        "bii",
        "broad_type",
        "rmag",
        "bmag",
        "optical_flag",
        "red_psf_flag",
        "blue_psf_flag",
        "redshift",
        "ref_name",
        "ref_redshift",
        "qso_prob",
        "radio_name",
        "xray_name",
        "alt_name_1",
        "alt_name_2",
        "class",
    ]
    mqp.assign_file_reader(
        reader_func=pd.read_table,
        read_chunks=True,
        names=columns,
        chunksize=50000,
        engine="c",
        skiprows=65,
        skip_blank_lines=True,
        sep="|",
        index_col=False,
        comment="<",
    )

    def mqc_modifier(srcdict):
        # format coordinates into geoJSON type (commented version uses 'legacy' pair):
        # unfortunately mongo needs the RA to be folded into -180, +180
        ra = srcdict["ra"] if srcdict["ra"] < 180.0 else srcdict["ra"] - 360.0
        srcdict["pos"] = {"type": "Point", "coordinates": [ra, srcdict["dec"]]}

        # add healpix index
        srcdict["hpxid_16"] = int(
            ang2pix(2 ** 16, srcdict["ra"], srcdict["dec"], lonlat=True, nest=True)
        )

        return srcdict

    mqp.assign_dict_modifier(mqc_modifier)

    mqp.push_to_db(
        coll_name="srcs",
        index_on=["hpxid_16", [("pos", pymongo.GEOSPHERE)]],
        index_args=[{}, {}],  # specify arguments for the index creation
        overwrite_coll=False,
        append_to_coll=False,
        dbclient=mongo_client,
    )

    assert mqp.coll.estimated_document_count() == 34

    mqp.healpix_meta(healpix_id_key="hpxid_16", order=16, is_indexed=True, nest=True)
    mqp.sphere2d_meta(sphere2d_key="pos", is_indexed=True, pos_format="geoJSON")
    mqp.science_meta(
        contact="C. Norris",
        email="chuck.norris@desy.de",
        description="compilation of AGN and Quasar",
        reference="http://quasars.org/milliquas.htm",
    )
