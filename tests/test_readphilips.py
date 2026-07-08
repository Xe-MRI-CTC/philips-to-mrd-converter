import readphilips.ReadPhilips as rp
import datasets
import pytest


@pytest.fixture(autouse=True)
def test_data():
    test_data = datasets.getTestData()
    return test_data


def test_dl_instantiate(test_data):
    for dataset in range(len(test_data)):
        dl_name = test_data[dataset]["dl"]
        result = rp.PhilipsData(dl_name)
        assert isinstance(result, rp.PhilipsData)


def test_rls_instantiate(test_data):
    for dataset in range(len(test_data)):
        rls_name = test_data[dataset]["rls"]
        result = rp.PhilipsData(rls_name)
        assert isinstance(result, rp.PhilipsData)


def test_dl_compute(test_data):
    for dataset in range(len(test_data)):
        dl_name = test_data[dataset]["dl"]
        result = rp.PhilipsData(dl_name)
        result.compute()
        assert isinstance(result, rp.PhilipsData)


def test_rls_compute(test_data):
    for dataset in range(len(test_data)):
        rls_name = test_data[dataset]["rls"]
        result = rp.PhilipsData(rls_name)
        assert isinstance(result, rp.PhilipsData)


if __name__ == "__main__":
    pytest.main()
