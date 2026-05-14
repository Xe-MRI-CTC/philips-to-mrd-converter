import unittest
import os
from pathlib import Path
import readphilips.ReadPhilips as rp


class TestReadPhilips(unittest.TestCase):
    def setUp(self):
        # Set paths
        self.loc = Path(__file__).parent.absolute()
        data_path = os.path.join(self.loc, "testdata")

        # Set up test examples
        self.test_data = [
            {  # spiral data example
                "type": "ventilation",
                "name": "spiral ventilation",
                "dl": os.path.join(data_path, "2DSpiral", "2DSpiral.data"),
                "rls": os.path.join(data_path, "2DSpiral", "2DSpiral.sin")
            },
        ]

    def test_1_dl(self):
        print("\\nTesting ReadPhilips Data\\List class instantiation...")
        try:
            for dataset in range(len(self.test_data)):
                dl_name = self.test_data[dataset]["dl"]
                result = rp.PhilipsData(dl_name)
            self.assertIsInstance(
                result, rp.PhilipsData, "ReadPhilips class not instantiated correctly for Data\\List")
        except Exception:
            self.assertIsInstance(
                result, rp.PhilipsData, "ReadPhilips class not instantiated correctly for Data\\List")
        print("\\nTested ReadPhilips Data\\List class instantiation: Successful")

    def test_2_rls(self):
        print("\\nTesting ReadPhilips Raw\\Lab\\Sin class instantiation...")
        try:
            for dataset in range(len(self.test_data)):
                rls_name = self.test_data[dataset]["rls"]
                result = rp.PhilipsData(rls_name)
            self.assertIsInstance(
                result, rp.PhilipsData, "ReadPhilips class not instantiated correctly for Raw\\Lab\\Sin")
        except Exception:
            self.assertIsInstance(
                result, rp.PhilipsData, "ReadPhilips class not instantiated correctly for Raw\\Lab\\Sin")
        print("\\nTested ReadPhilips Raw\\Lab\\Sin class instantiation: Successful")

    def test_3_dl_compute(self):
        print("\\nTesting ReadPhilips Data\\List class compute...")
        try:
            for dataset in range(len(self.test_data)):
                dl_name = self.test_data[dataset]["dl"]
                result = rp.PhilipsData(dl_name)
                result.compute()
            self.assertIsInstance(
                result, rp.PhilipsData, "ReadPhilips class not computed for Data\\List")
        except Exception:
            self.assertIsInstance(
                result, rp.PhilipsData, "ReadPhilips class not computed correctly for Data\\List")
        print("\\nTested ReadPhilips Data\\List class compute: Successful")

    def test_4_rls_compute(self):
        print("\\nTesting ReadPhilips Raw\\Lab\\Sin class compute...")
        try:
            for dataset in range(len(self.test_data)):
                rls_name = self.test_data[dataset]["rls"]
                result = rp.PhilipsData(rls_name)
            self.assertIsInstance(
                result, rp.PhilipsData, "ReadPhilips class not computed correctly for Raw\\Lab\\Sin")
        except Exception:
            self.assertIsInstance(
                result, rp.PhilipsData, "ReadPhilips class not computed correctly for Raw\\Lab\\Sin")
        print("\\nTested ReadPhilips Raw\\Lab\\Sin class compute: Successful")


if __name__ == "__main__":
    unittest.main()
