#%% Import modules:
import os

import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits import mplot3d

import readphilips as rp

# Working directory:
location = os.getcwd()
    
#%% 2D spiral data:
    
# Read in the raw-list files:
filename = location + '\\data\\2D Spiral Ventilation CCHMC\\' + 'raw_004.data' 

# Read in the data using the ReadPhilips.py script:
inputfile = rp.PhilipsData(filename)
inputfile.compute()

# Extract data:
data = abs(inputfile.data)

# Plot an FID to check data:
fid = data[0,1,1,:] # First dynamic, first slice, first spiral
plt.figure()
plt.plot(fid)
plt.ylabel('k-space intensity')
plt.xlabel('Sample number')
plt.title('Free induction decay')
plt.show()

# Read in the raw-lab-sin files:
filename = location + '\\data\\2D Spiral Ventilation CCHMC\\' + '20211013_113717_CPIR_Vent_HANNING_2DSOS_WIP.sin' 

# Read in the data using the ReadPhilips.py script:
inputfile =  rp.PhilipsData(filename)
inputfile.compute()

# Extract spiral coordinates:
coords = inputfile.spparams.get('COORDS')

# Plot the coordinates:
plt.figure()
plt.plot(coords[0,:], coords[1,:])
plt.ylabel('$k_y$')
plt.xlabel('$k_x$')
plt.title('Spiral coordinates')
plt.show()

# Extract expanded spiral coordinates:
coords_expanded = inputfile.spparams.get('COORDS_EXPANDED')

# Plot the coordinates:
plt.figure()
plt.plot(coords_expanded[0,:,:], coords_expanded[1,:,:])
plt.ylabel('$k_y$')
plt.xlabel('$k_x$')
plt.title('Spiral coordinates')
plt.show()

#%% 3D radial data:
    
# Read in the raw-list files:
filename = location + '\\data\\3D Radial Gas-exchange CCHMC\\' + 'raw_007.list' 
# Note: the Duke data is in the correct radial format, whereas the CCHMC data 
# has had some editing so the data comes out strange.

# For future reference, CCHMC data looks like this:
    # "Approximate translation" from complicated Philips settings to MRI:
    # Dim 1 = mixes; used to separate spectroscopy data and imaging data.
    # Dim 2 = dynamics; used to separate gas, dissolved, and off-res acquisitions.
    # Dim 3 = kz; Archimedean spiral interleaves.
    # Dim 4 = ky; projections per interleave.
    # Dim 5 = kx; read out.
    # Useful paper: https://onlinelibrary.wiley.com/doi/full/10.1002/mrm.22898

# Read in the data using the ReadPhilips.py script:
inputfile =  rp.PhilipsData(filename)
inputfile.compute()

# Extract data:
data = inputfile.data

# Plot an FID to check data:
fid = abs(data[0,0,0,0,0:58])

plt.figure()
plt.plot(fid)
plt.ylabel('k-space intensity')
plt.xlabel('Sample number')
plt.title('Free induction decay')
plt.show()

# Read in the raw-lab-sin files:
filename = location + '\\data\\3D Radial Gas-exchange CCHMC\\' + '20191008_162511_Dissolved_Xe_20191008.sin'

# Read in the data using the ReadPhilips.py script:
inputfile =  rp.PhilipsData(filename)
inputfile.readParamOnly = True # Read the sin file only, don't bother with raw-lab. 
inputfile.compute()

# Extract radial coordinates:
coords = inputfile.radparams.get('COORDS')

# Plot the coordinates:
plt.figure()
ax = plt.axes(projection='3d')
ax.scatter(coords[:,:,-1,0], coords[:,:,-1,1], coords[:,:,-1,2], c = 'r', s = 1)
ax.set_zlabel('$k_z$')
ax.set_ylabel('$k_y$')
ax.set_xlabel('$k_x$')
plt.title('Radial coordinates')
plt.show()


# Extract data for recon:
data1 = data[0,0,:,:,:58]
data2 = np.transpose(data1.reshape(data1.shape[-1], -1))
data = data2
coords1 = np.reshape(coords, (950,58,3))
traj = coords1
