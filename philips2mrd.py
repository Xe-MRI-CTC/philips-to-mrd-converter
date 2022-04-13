import os
import ismrmrd as mrd
from rp.readphilips import *

loc = os.getcwd()
fname = loc + '\\rp\\data\\raw_001.data' #2dim

out = mrd.Acquisition()

out = mrd.Acquisition()

inputfile = PhilipsData(fname)
inputfile.compute()

