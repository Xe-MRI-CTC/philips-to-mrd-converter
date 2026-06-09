import unittest
from unittest.mock import MagicMock
import sys
from pathlib import Path
import os

# Mock tkinter for github actions runners
sys.modules['tkinter'] = MagicMock()
sys.modules['tkinter.filedialog'] = MagicMock()

from philips2mrd import Ph2Mrd  # noqa: E402
from Scripts import XeGasExchange2XeCTCMRD  # noqa: E402


class TestPhilips2MRD(unittest.TestCase):
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
                "rls": os.path.join(data_path, "2DSpiral", "2DSpiral.sin"),
                "traj": None
            },
            {  # ctc gas ex example
                "type": "gas exchange",
                "name": "CTC gas exchange",
                "dl": os.path.join(data_path, "3DRadial_GXCTC", "3DRadial_GXCTC.data"),
                "rls": os.path.join(data_path, "3DRadial_GXCTC", "3DRadial_GXCTC.sin"),
                "traj": None
            },
            {  # gas ex w/ bonus spectra example
                "type": "gas exchange",
                "name": "gas exchange w/ bonus spectra",
                "dl": os.path.join(data_path, "3DRadial_GXBonus", "3DRadial_GXBonus.data"),
                "rls": os.path.join(data_path, "3DRadial_GXBonus", "3DRadial_GXBonus.sin"),
                "traj": os.path.join(data_path, "3DRadial_GXBonus", "3DRadial_GXBonus_Traj.sin")
            },
            {  # gas ex w/ 2 echo example
                "type": "gas exchange",
                "name": "gas exchange w/ 2 echoes",
                "dl": os.path.join(data_path, "3DRadial_GX2Echo", "3DRadial_GX2Echo.data"),
                "rls": os.path.join(data_path, "3DRadial_GX2Echo", "3DRadial_GX2Echo.sin"),
                "traj": None
            },
            {  # gas ex floret example
                "type": "gas exchange",
                "name": "gas exchange w/ floret",
                "dl": os.path.join(data_path, "3DFLORET_GXBonus", "3DFLORET_GXBonus.data"),
                "rls": os.path.join(data_path, "3DFLORET_GXBonus", "3DFLORET_GXBonus.sin"),
                "traj": os.path.join(data_path, "3DFLORET_GXBonus", "3DFLORET_GXBonus_Traj.sin")
            },
            {  # XeCTC Calibration example
                "type": "calibration",
                "name": "CTC calibration",
                "dl": os.path.join(data_path, "XeCTC_Calibration", "XeCTC_Calibration.data"),
                "rls": os.path.join(data_path, "XeCTC_Calibration", "XeCTC_Calibration.sin"),
                "traj": None
            },
        ]

    def test_1_setup(self):
        print("\\nTesting empty Ph2Mrd class instantiation...")
        result = Ph2Mrd()
        self.assertIsInstance(
            result, Ph2Mrd, "Empty Ph2Mrd Class not instantiated correctly")

    def test_2_both(self):
        print("\\nTesting class instantiation with both raw data files...")
        for dataset in range(len(self.test_data)):
            data_name = self.test_data[dataset]["name"]
            # create converter instance
            print(
                f"    Testing Ph2Mrd class instantiation with both raw types for data set {dataset} : {data_name}")
            result = Ph2Mrd(self.test_data[dataset]
                            ["dl"], self.test_data[dataset]["rls"])
            # confirm instance created
            self.assertIsInstance(
                result, Ph2Mrd, "Full Ph2Mrd Class not instantiated correctly")
            # run converter
            print(f"    Testing data set {dataset} conversion")
            result.convert(self.loc)

    def test_3_dl(self):
        print("\\nTesting class instantiation with only data/list files...")
        for dataset in range(len(self.test_data)):
            data_name = self.test_data[dataset]["name"]
            print(
                f"    Testing Ph2Mrd class instantiation with only data/list for data set {dataset} : {data_name}")
            # create converter instance
            result = Ph2Mrd(self.test_data[dataset]["dl"], None)
            # confirm instance created
            self.assertIsInstance(
                result, Ph2Mrd, "Full Ph2Mrd Class not instantiated correctly")
            # run converter
            print(f"    Testing data set {dataset} conversion")
            result.convert(self.loc)

    def test_4_rls(self):
        print("\\nTesting class instantiation with only raw/lab/sin files...")
        for dataset in range(len(self.test_data)):
            data_name = self.test_data[dataset]["name"]
            print(
                f"    Testing Ph2Mrd class instantiation with only raw/lab/sin for data set {dataset} : {data_name}")
            # create converter instance
            result = Ph2Mrd(None, self.test_data[dataset]["rls"])
            # confirm instance created
            self.assertIsInstance(
                result, Ph2Mrd, "Full Ph2Mrd Class not instantiated correctly")
            # run converter
            print(f"    Testing data set:{dataset} conversion")
            result.convert(self.loc)

    def test_5_gx_converter_both(self):
        print("\\nTesting gx converter script with both raw data files...")
        for dataset in range(len(self.test_data)):
            data_name = self.test_data[dataset]["name"]
            if self.test_data[dataset]["type"] == "gas exchange":
                print(
                    f"    Testing gx converter script for both raw data for data set {dataset} : {data_name}")
                XeGasExchange2XeCTCMRD.Gx2XeCTCMRD(
                    data_file=self.test_data[dataset]["dl"], raw_file=self.test_data[dataset]["rls"],
                    traj_file=self.test_data[dataset]["traj"])

    def test_6_gx_converter_rls(self):
        print("\\nTesting gx converter script with only rls data files...")
        for dataset in range(len(self.test_data)):
            data_name = self.test_data[dataset]["name"]
            if self.test_data[dataset]["type"] == "gas exchange":
                print(
                    f"    Testing gx converter script for only rls data for data set {dataset} : {data_name}")
                XeGasExchange2XeCTCMRD.Gx2XeCTCMRD(
                    data_file='', raw_file=self.test_data[dataset]["rls"], traj_file=self.test_data[dataset]["traj"])

    def test_7_gx_converter_dl(self):
        print("\\nTesting gx converter script with only dl data files...")
        for dataset in range(len(self.test_data)):
            data_name = self.test_data[dataset]["name"]
            if self.test_data[dataset]["type"] == "gas exchange":
                print(
                    f"    Testing gx converter script for only dl data for data set {dataset} : {data_name}")
                XeGasExchange2XeCTCMRD.Gx2XeCTCMRD(
                    data_file=self.test_data[dataset]["dl"], raw_file='', traj_file=self.test_data[dataset]["traj"])

    def test_8_cal_converter_both(self):
        print("\\nTesting gx calibration converter script with both raw data files...")
        for dataset in range(len(self.test_data)):
            data_name = self.test_data[dataset]["name"]
            if self.test_data[dataset]["type"] == "calibration":
                print(
                    f"    Testing gx calibration converter script for both raw data for data set {dataset} : {data_name}")
                XeGasExchange2XeCTCMRD.Gx2XeCTCMRD(
                    data_file=self.test_data[dataset]["dl"], raw_file=self.test_data[dataset]["rls"])


if __name__ == "__main__":
    unittest.main()
