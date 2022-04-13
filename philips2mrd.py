import os
from rp.ReadPhilips import PhilipsData

loc = os.getcwd()
fname = loc + '\\rp\\data\\raw_001.data' #2dim

<<<<<<< Updated upstream
=======
#fname = 'fakeloc_asdaskdfjak'
os.path.exists(fname)
out = mrd.Acquisition()

>>>>>>> Stashed changes
inputfile = PhilipsData(fname)
inputfile.compute()

