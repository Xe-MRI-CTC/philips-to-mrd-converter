import pytest
from pathlib import Path
import datasets
from philips2mrd.philips2mrd import Ph2Mrd
from Scripts import XeGasExchange2XeCTCMRD


@pytest.fixture(autouse=True)
def test_data():
    test_data = datasets.getTestData()
    return test_data


@pytest.fixture(autouse=True)
def location():
    loc = Path(__file__).parent.absolute()
    return loc


@pytest.fixture(autouse=True)
def cleanup(location):
    yield
    for file_path in location.rglob("*.h5"):
        if file_path.is_file():  # Safety check to ensure it's a file
            print('\nDeleting ' + str(file_path))
            file_path.unlink()


def test_instantiate_empty():
    result = Ph2Mrd()
    assert isinstance(result, Ph2Mrd), "Empty Ph2Mrd instantiation failed..."


def test_convert_both(test_data, location, subtests):
    for dataset in range(len(test_data)):
        data_name = test_data[dataset]["name"]
        # create converter instance
        result = Ph2Mrd(test_data[dataset]["dl"], test_data[dataset]["rls"])
        # confirm instance created
        with subtests.test(msg="[Both] Testing instantiation: " + data_name, i=dataset):
            assert isinstance(result, Ph2Mrd), "[Both] Ph2Mrd instantiation failed..."
        # run converter
        mrdFileName, rlsPhData, dlPhData = result.convert(location)
        # confirm instance converted
        with subtests.test(msg="[Both] Testing conversion: " + data_name, i=dataset):
            assert mrdFileName.is_file(), "[Both] Ph2Mrd conversion failed..."
        with subtests.test(msg="[Both] Checking rls return: " + data_name, i=dataset):
            assert rlsPhData is not None, "[Both] rls data return failed..."
        with subtests.test(msg="[Both] Checking dl return: " + data_name, i=dataset):
            assert dlPhData is not None, "[Both] dl data return failed..."


def test_convert_rls(test_data, location, subtests):
    for dataset in range(len(test_data)):
        data_name = test_data[dataset]["name"]
        # create converter instance
        result = Ph2Mrd(None, test_data[dataset]["rls"])
        # confirm instance created
        with subtests.test(msg="[RLS] Testing instantiation: " + data_name, i=dataset):
            assert isinstance(result, Ph2Mrd), "[RLS] Ph2Mrd instantiation failed..."
        # run converter
        mrdFileName, rlsPhData, dlPhData = result.convert(location)
        # confirm instance converted
        with subtests.test(msg="[RLS] Testing conversion: " + data_name, i=dataset):
            assert mrdFileName.is_file(), "[RLS] Ph2Mrd conversion failed..."
        with subtests.test(msg="[RLS] Checking rls return: " + data_name, i=dataset):
            assert rlsPhData is not None, "[RLS] rls data return failed..."
        with subtests.test(msg="[RLS] Checking no dl return: " + data_name, i=dataset):
            assert dlPhData is None, "[RLS] dl data return unexpected..."


def test_convert_dl(test_data, location, subtests):
    for dataset in range(len(test_data)):
        data_name = test_data[dataset]["name"]
        # create converter instance
        result = Ph2Mrd(test_data[dataset]["dl"], None)
        # confirm instance created
        with subtests.test(msg="[DL] Testing instantiation: " + data_name, i=dataset):
            assert isinstance(result, Ph2Mrd), "[DL] Ph2Mrd instantiation failed..."
        # run converter
        mrdFileName, rlsPhData, dlPhData = result.convert(location)
        # confirm instance converted
        with subtests.test(msg="[DL] Testing conversion: " + data_name, i=dataset):
            assert mrdFileName.is_file(), "[DL] Ph2Mrd conversion failed..."
        with subtests.test(msg="[DL] Checking no rls return: " + data_name, i=dataset):
            assert rlsPhData is None, "[DL] rls data return unexpected..."
        with subtests.test(msg="[DL] Checking dl return: " + data_name, i=dataset):
            assert dlPhData is not None, "[RLS] dl data return failed..."


def test_gx_converter_both(test_data, subtests):
    for dataset in range(len(test_data)):
        data_name = test_data[dataset]["name"]
        if test_data[dataset]["type"] != "standard":
            output_path = XeGasExchange2XeCTCMRD.Gx2XeCTCMRD(
                data_file=test_data[dataset]["dl"], raw_file=test_data[dataset]["rls"],
                traj_file=test_data[dataset]["traj"])
            with subtests.test(msg="[GX Both] Testing conversion: " + data_name, i=dataset):
                assert output_path.is_file(), "[GX Both] Ph2Mrd conversion failed..."


def test_gx_converter_rls(test_data, subtests):
    for dataset in range(len(test_data)):
        data_name = test_data[dataset]["name"]
        if test_data[dataset]["type"] != "standard":
            output_path = XeGasExchange2XeCTCMRD.Gx2XeCTCMRD(
                raw_file=test_data[dataset]["rls"],
                traj_file=test_data[dataset]["traj"])
            with subtests.test(msg="[GX RLS] Testing conversion: " + data_name, i=dataset):
                assert output_path.is_file(), "[GX RLS] Ph2Mrd conversion failed..."


def test_gx_converter_dl(test_data, subtests):
    for dataset in range(len(test_data)):
        data_name = test_data[dataset]["name"]
        if test_data[dataset]["type"] != "standard":
            output_path = XeGasExchange2XeCTCMRD.Gx2XeCTCMRD(
                data_file=test_data[dataset]["dl"],
                traj_file=test_data[dataset]["traj"])
            with subtests.test(msg="[GX DL] Testing conversion: " + data_name, i=dataset):
                assert output_path.is_file(), "[GX DL] Ph2Mrd conversion failed..."


if __name__ == "__main__":
    pytest.main()
