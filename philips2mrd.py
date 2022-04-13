import os
from rp.ReadPhilips import PhilipsData

loc = os.getcwd()
fname = loc + '\\rp\\data\\raw_001.data' #2dim

out = mrd.Acquisition()

inputfile = PhilipsData(fname)
inputfile.compute()

