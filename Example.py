from pathlib import Path
import numpy as np
import ismrmrd as mrd
import ismrmrd.xsd
import philips2mrd as p2m
import matplotlib.pyplot as plt

# Set paths
loc = Path(__file__).parent.absolute()
dlName = loc / "testdata" / "2DSpiral" / "2DSpiral.data"
rlsName = loc / "testdata" / "2DSpiral" / "2DSpiral.sin"
outDir = loc

# Run converter
inputData = p2m.Ph2Mrd(dlName, rlsName)
inputData.trajorder = 2
inputData.delay = -1.25
mrdName, rls, dl = inputData.convert(outDir)

# Load
dset = mrd.Dataset(mrdName, "dataset", create_if_needed=False)

# Get info
header = ismrmrd.xsd.CreateFromDocument(dset.read_xml_header())
enc = header.encoding[0]

# Matrix size
eNx = enc.encodedSpace.matrixSize.x
eNy = enc.encodedSpace.matrixSize.y
eNz = enc.encodedSpace.matrixSize.z

ncoils = header.acquisitionSystemInformation.receiverChannels
if enc.encodingLimits.slice != None:
    nslices = enc.encodingLimits.slice.maximum + 1
else:
    nslices = 1

if enc.encodingLimits.repetition != None:
    nreps = enc.encodingLimits.repetition.maximum + 1
else:
    nreps = 1

if enc.encodingLimits.contrast != None:
    ncontrasts = enc.encodingLimits.contrast.maximum + 1
else:
    ncontrasts = 1

if enc.encodingLimits.kspace_encoding_step_2.maximum == 0:
    ndim = 2
else:
    ndim = 3

# loop through the acquisitions looking for noise scans
firstacq = 0
for acqnum in range(dset.number_of_acquisitions()):
    acq = dset.read_acquisition(acqnum)

    # Currently ignoring noise scans
    if acq.isFlagSet(ismrmrd.ACQ_IS_NOISE_MEASUREMENT):
        print("Found noise scan at acq ", acqnum)
        continue
    else:
        firstacq = acqnum
        print("Imaging acquisition starts acq ", acqnum)
        break

# Initialiaze a storage array
all_data = np.zeros((nreps, ncontrasts, nslices, ncoils,
                    eNz, eNy, eNx), dtype=np.complex64)
all_traj = np.zeros((nreps, ncontrasts, nslices,
                    eNz, eNy, eNx, ndim), dtype=np.float32)

# Loop through the rest of the acquisitions and stuff
for acqnum in range(firstacq, dset.number_of_acquisitions()):
    acq = dset.read_acquisition(acqnum)

    # Stuff into the buffer
    rep = acq.idx.repetition
    contrast = acq.idx.contrast
    slice = acq.idx.slice
    y = acq.idx.kspace_encode_step_1
    z = acq.idx.kspace_encode_step_2
    all_data[rep, contrast, slice, :, z, y, :] = acq.data
    all_traj[rep, contrast, slice, z, y, :, :] = acq.traj

print(dset.read_xml_header().decode('utf-8'))

dset.close()

# Plot trajectory
fig = plt.figure(figsize=(13, 6))
traj = all_traj[0,0,0,0,:,:,:]
n_traj = 100
n_traj = min(n_traj, len(traj))
if acq.trajectory_dimensions == 2:
    ax = fig.add_subplot(121)
else:
    ax = fig.add_subplot(121, projection="3d")
colors = plt.cm.jet(np.linspace(0,1,n_traj))
for i in range(n_traj):
    ax.plot(*traj[i].T, color=colors[i], lw=1)
ax.set_aspect('equal')
ax.set_title('Trajectory')

# Plot data
slice_plot = nslices//2
data = all_data[0,0,slice_plot,0,0,:,:]
n_data = 100
n_data = min(n_data, len(data))
ax = fig.add_subplot(122)
colors = plt.cm.jet(np.linspace(0,1,n_traj))
for i in range(n_data):
    ax.plot(np.abs(data[i,:]), color=colors[i], lw=1)
ax.set_title('Data')
plt.show()

print("Data converted to mrd format: ", mrdName)
