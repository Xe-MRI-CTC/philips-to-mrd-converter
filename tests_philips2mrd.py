import unittest
import os
from philips2mrd import *
from Scripts import XeGasExchange2XeCTCMRD


class TestPhilips2MRD(unittest.TestCase):
    def setUp(self):
        # Set paths
        self.loc = Path(__file__).parent.absolute()
        data_path = os.path.join(self.loc, "rp", "data")

        # Set up test examples
        self.test_data = [
            {  # spiral data example
                "type": "ventilation",
                "name": "spiral ventilation",
                "dl": os.path.join(data_path, "2DSpiralVentilationCCHMC", "raw_004.data"),
                "rls": os.path.join(data_path, "2DSpiralVentilationCCHMC", "20211013_113717_CPIR_Vent_HANNING_2DSOS_WIP.sin")
            },
            {  # ctc gas ex example
                "type": "gas exchange",
                "name": "CTC gas exchange",
                "dl": os.path.join(data_path, "3DRadialGas-exchangeDuke", "raw_207.data"),
                "rls": os.path.join(data_path, "3DRadialGas-exchangeDuke", "20220112_113653_DukeIPF_Gas_Exchange.sin")
            },
            {  # gas ex w/ bonus spectra example
                "type": "gas exchange",
                "name": "gas exchange w/ bonus spectra",
                "dl": os.path.join(data_path, "3DRadialGas-exchangeCCHMC", "raw_007.data"),
                "rls": os.path.join(data_path, "3DRadialGas-exchangeCCHMC", "20191008_162511_Dissolved_Xe_20191008.sin")
            },
            {  # gas ex w/ 2 echo example
                "type": "gas exchange",
                "name": "gas exchange w/ 2 echoes",
                "dl": os.path.join(data_path, "3DRadialGas-exchange2Echo", "raw_1501.data"),
                "rls": os.path.join(data_path, "3DRadialGas-exchange2Echo", "20250218_134601_Xenon_3D_radial_1ptDixon_2Echoes.sin")
            },
            {  # gas ex floret example
                "type": "gas exchange",
                "name": "gas exchange w/ floret",
                "dl": os.path.join(data_path, "3DFLORETGas-exchange", "raw_111.data"),
                "rls": os.path.join(data_path, "3DFLORETGas-exchange", "20251002_152616_Xenon_3D_FLORET_Dixon.sin")
            },
        ]

    def test_setup(self):
        print(f"Testing empty Ph2Mrd class instantiation...")
        result = Ph2Mrd()
        self.assertIsInstance(
            result, Ph2Mrd, "Empty Ph2Mrd Class not instantiated correctly")

    def test_both(self):
        print(f"Testing class instantiation with both raw data files...")
        for dataset in range(len(self.test_data)):
            data_name = self.test_data[dataset]["name"]
            # create converter instance
            print(f"    Testing Ph2Mrd class instantiation with both raw types for data set {dataset} : {data_name}")
            result = Ph2Mrd(self.test_data[dataset]
                            ["dl"], self.test_data[dataset]["rls"])
            # confirm instance created
            self.assertIsInstance(
                result, Ph2Mrd, "Full Ph2Mrd Class not instantiated correctly")
            # run converter
            print(f"    Testing data set {dataset} conversion")
            result.convert(self.loc)

    def test_dl(self):
        print(f"Testing class instantiation with only data/list files...")
        for dataset in range(len(self.test_data)):
            data_name = self.test_data[dataset]["name"]
            print(f"    Testing Ph2Mrd class instantiation with only data/list for data set {dataset} : {data_name}")
            # create converter instance
            result = Ph2Mrd(self.test_data[dataset]["dl"], None)
            # confirm instance created
            self.assertIsInstance(
                result, Ph2Mrd, "Full Ph2Mrd Class not instantiated correctly")
            # run converter
            # print(f"    Testing data set {dataset} conversion")
            # result.convert(self.loc)

    def test_rls(self):
        print(f"Testing class instantiation with only raw/lab/sin files...")
        for dataset in range(len(self.test_data)):
            data_name = self.test_data[dataset]["name"]
            print(f"    Testing Ph2Mrd class instantiation with only raw/lab/sin for data set {dataset} : {data_name}")
            # create converter instance
            result = Ph2Mrd(None, self.test_data[dataset]["rls"])
            # confirm instance created
            self.assertIsInstance(
                result, Ph2Mrd, "Full Ph2Mrd Class not instantiated correctly")
            # run converter
            # print(f"    Testing data set:{dataset} conversion")
            # result.convert(self.loc)

    def test_gx_converter_script(self):
        print(f"Testing gx converter script with both raw data files...")
        for dataset in range(len(self.test_data)):
            data_name = self.test_data[dataset]["name"]
            if self.test_data[dataset]["type"] == "gas exchange":
                print(f"    Testing gx converter script for both raw data for data set {dataset} : {data_name}")
                XeGasExchange2XeCTCMRD.Gx2XeCTCMRD(
                    data_file=self.test_data[1]["dl"], raw_file=self.test_data[1]["rls"], traj_file=None)


if __name__ == "__main__":
    unittest.main()
