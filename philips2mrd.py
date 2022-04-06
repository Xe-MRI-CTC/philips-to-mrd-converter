import os
from rp.ReadPhilips import PhilipsData

loc = os.getcwd()
fname = loc + '\\data\\raw_001.data' #2dim

inputfile = PhilipsData(fname)
inputfile.compute()

