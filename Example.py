
import os
from ReadPhilips import PhilipsData
from matplotlib import pyplot as plt

loc = os.getcwd()
#fname = loc + '\\data\\raw_001.data' #2dim
fname = loc + '\\data\\20220112_113653_DukeIPF_Gas_Exchange.sin' #3dim

inputfile = PhilipsData(fname)
inputfile.compute()

#dataSlice = abs(inputfile.data)
dataSlice = abs(inputfile.data[1,1:64,:])

im = plt.imshow(dataSlice, cmap='hot')
plt.show()
