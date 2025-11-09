import unittest
from philips2mrd import *


class TestPhilips2MRD(unittest.TestCase):
    def setUp(self):
        # Set paths
        self.loc = Path(__file__).parent.absolute()
        self.dl_gx = self.loc / "rp" / "data" / \
            "3DRadialGas-exchangeDuke" / "raw_207.data"
        self.rls_gx = self.loc / "rp" / "data" / "3DRadialGas-exchangeDuke" / \
            "20220112_113653_DukeIPF_Gas_Exchange.sin"
        self.dl_gx_bonus = self.loc / "rp" / "data" / \
            "3DRadialGas-exchangeCCHMC" / "raw_007.data"
        self.rls_gx_bonus = self.loc / "rp" / "data" / "3DRadialGas-exchangeCCHMC" / \
            "20191008_162511_Dissolved_Xe_20191008.sin"

    def test_setup(self):
        result = Ph2Mrd()
        self.assertIsInstance(
            result, Ph2Mrd, "Empty Ph2Mrd Class not instantiated correctly")

    def test_gx_both(self):
        result = Ph2Mrd(self.dl_gx, self.rls_gx)
        self.assertIsInstance(
            result, Ph2Mrd, "Full Ph2Mrd Class not instantiated correctly")
        result.convert(self.loc)

    def test_gx_dl(self):
        result = Ph2Mrd(self.dl_gx, None)
        self.assertIsInstance(
            result, Ph2Mrd, "Data\list Ph2Mrd Class not instantiated correctly")
        # result.convert(self.loc)

    def test_gx_rls(self):
        result = Ph2Mrd(None, self.rls_gx)
        self.assertIsInstance(
            result, Ph2Mrd, "Raw\Lab\Sin Ph2Mrd Class not instantiated correctly")
        # result.convert(self.loc)

    def test_gx_bonus_both(self):
        result = Ph2Mrd(self.dl_gx_bonus, self.rls_gx_bonus)
        self.assertIsInstance(
            result, Ph2Mrd, "Full Ph2Mrd Class not instantiated correctly")
        result.convert(self.loc)

    def test_gx_bonus_dl(self):
        result = Ph2Mrd(self.dl_gx_bonus, None)
        self.assertIsInstance(
            result, Ph2Mrd, "Data\list Ph2Mrd Class not instantiated correctly")
        # result.convert(self.loc)

    def test_gx_bonus_rls(self):
        result = Ph2Mrd(None, self.rls_gx_bonus)
        self.assertIsInstance(
            result, Ph2Mrd, "Raw\Lab\Sin Ph2Mrd Class not instantiated correctly")
        # result.convert(self.loc)


if __name__ == '__main__':
    unittest.main()
