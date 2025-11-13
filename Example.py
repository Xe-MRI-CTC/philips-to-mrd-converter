from pathlib import Path
import numpy as np
import ismrmrd as mrd
import ismrmrd.xsd
import philips2mrd as p2m

# Set paths
loc = Path(__file__).parent.absolute()
dlName = loc / "rp" / "data" / "3DFLORETGas-exchange" / "raw_111.data"
rlsName = loc / "rp" / "data" / "3DFLORETGas-exchange" / \
    "20251002_152616_Xenon_3D_FLORET_Dixon.sin"
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

# Loop through the rest of the acquisitions and stuff
for acqnum in range(firstacq, dset.number_of_acquisitions()):
    acq = dset.read_acquisition(acqnum)

    if acq.idx.contrast > 0 or acq.idx.set > 0: # Skip those which may have different readout lengths
        continue

    # Stuff into the buffer
    rep = acq.idx.repetition
    contrast = acq.idx.contrast
    slice = acq.idx.slice
    y = acq.idx.kspace_encode_step_1
    z = acq.idx.kspace_encode_step_2
    all_data[rep, contrast, slice, :, z, y, :] = acq.data

dset.close()

print("Data converted to mrd format: ", mrdName)
