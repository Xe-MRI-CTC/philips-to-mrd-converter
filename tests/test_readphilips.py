import readphilips.ReadPhilips as rp
import datasets
import pytest
from philips2mrd.functions import mrd_recon as mrr


@pytest.fixture(autouse=True)
def test_data():
    test_data = datasets.getTestData()
    return test_data


def test_rls(test_data, subtests):
    for dataset in range(len(test_data)):
        rls_name = test_data[dataset]["rls"]
        result = rp.PhilipsData(rls_name)
        # confirm instantiation
        with subtests.test(msg="[RLS] Testing instantiation: " + rls_name, i=dataset):
            assert isinstance(result, rp.PhilipsData), "[RLS] RP instantiation failed..."
        # confirm compute
        result.compute()
        # assert to confirm result.dat exists
        with subtests.test(msg="[RLS] Testing compute data: " + rls_name, i=dataset):
            assert hasattr(result, "data"), "[RLS] Data does not exist"
        # assert to confirm result.header exists
        with subtests.test(msg="[RLS] Testing compute header: " + rls_name, i=dataset):
            assert hasattr(result, "header"), "[RLS] Header does not exist"
        # confirm recon
        try:
            traj_type = int(result.header["sin"]["k_space_traj_type"][0][0])
        except Exception:
            traj_type = 0
        if traj_type == 0 and test_data[dataset]["type"] == 'standard':
            image = mrr.recon_cart_rp(result)
        elif test_data[dataset]["type"] == 'standard':
            image = mrr.recon_noncart_rp(result)
            continue
        else:
            continue
        with subtests.test(msg="[RLS] Testing image reconstruction: " + rls_name, i=dataset):
            assert 0.85 * test_data[dataset]["img_max"] <= image.max() <= 1.15 * test_data[dataset]["img_max"], "[RLS] Image does not match expectations"
            assert 0.85 * test_data[dataset]["img_min"] <= image.min() <= 1.15 * test_data[dataset]["img_min"], "[RLS] Image does not match expectations"


# TODO create DL test

# TODO create combo test


if __name__ == "__main__":
    pytest.main()
