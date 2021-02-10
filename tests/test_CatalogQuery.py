import pytest
from extcats.CatalogQuery import CatalogQuery


@pytest.fixture
def catq(milliquas, mongo_client):
    return CatalogQuery("milliquas", dbclient=mongo_client)


with_index_methods = pytest.mark.parametrize(
    "method",
    ["2dsphere", "healpix"],
)


@with_index_methods
def test_findwithin(catq, method):
    assert len(catq.findwithin(185.453, -89.241, 10, method=method)) == 1
    assert catq.findwithin(0, -90, 10, method=method) is None


@with_index_methods
def test_findclosest(catq, method):
    match, dist = catq.findclosest(185.453, -89.241, 10, method=method)
    assert dist == pytest.approx(3.5758, abs=1e-4)
    assert catq.findwithin(0, -90, 10) is None


@with_index_methods
def test_binarysearch(catq, method):
    assert catq.binaryserach(185.453, -89.241, 10, method=method) is True
    assert catq.binaryserach(0, -90, 10, method=method) is False
