import unittest
import os
from philips2mrd import *
from Scripts import XeGasExchange2XeCTCMRD


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
        ]

    def test_1_setup(self):
        print(f"\nTesting empty Ph2Mrd class instantiation...")
        result = Ph2Mrd()
        self.assertIsInstance(
            result, Ph2Mrd, "Empty Ph2Mrd Class not instantiated correctly")

    def test_2_both(self):
        print(f"\nTesting class instantiation with both raw data files...")
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
        print(f"\nTesting class instantiation with only data/list files...")
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
        print(f"\nTesting class instantiation with only raw/lab/sin files...")
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
        print(f"\nTesting gx converter script with both raw data files...")
        for dataset in range(len(self.test_data)):
            data_name = self.test_data[dataset]["name"]
            if self.test_data[dataset]["type"] == "gas exchange":
                print(
                    f"    Testing gx converter script for both raw data for data set {dataset} : {data_name}")
                XeGasExchange2XeCTCMRD.Gx2XeCTCMRD(
                    data_file=self.test_data[dataset]["dl"], raw_file=self.test_data[dataset]["rls"], traj_file=self.test_data[dataset]["traj"])

    def test_6_gx_converter_rls(self):
        print(f"\nTesting gx converter script with only rls data files...")
        for dataset in range(len(self.test_data)):
            data_name = self.test_data[dataset]["name"]
            if self.test_data[dataset]["type"] == "gas exchange":
                print(
                    f"    Testing gx converter script for only rls data for data set {dataset} : {data_name}")
                XeGasExchange2XeCTCMRD.Gx2XeCTCMRD(
                    data_file='', raw_file=self.test_data[dataset]["rls"], traj_file=self.test_data[dataset]["traj"])
                
    def test_7_gx_converter_dl(self):
        print(f"\nTesting gx converter script with only dl data files...")
        for dataset in range(len(self.test_data)):
            data_name = self.test_data[dataset]["name"]
            if self.test_data[dataset]["type"] == "gas exchange":
                print(
                    f"    Testing gx converter script for only dl data for data set {dataset} : {data_name}")
                XeGasExchange2XeCTCMRD.Gx2XeCTCMRD(
                    data_file=self.test_data[dataset]["dl"], raw_file='', traj_file=self.test_data[dataset]["traj"])                


if __name__ == "__main__":
    unittest.main()
