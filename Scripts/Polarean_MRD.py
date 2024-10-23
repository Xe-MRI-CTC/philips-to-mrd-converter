import sys
from pathlib import Path

p2mDir = Path(__file__).parent.parent.absolute()
sys.path.append(str(p2mDir))

import os
import datetime
import copy
import tkinter as tk
from tkinter import filedialog

import ismrmrd as mrd

import philips2mrd as p2m

# Get paths
root = tk.Tk()
root.withdraw()
dlName = filedialog.askopenfilename(title='Select .data file',filetypes=[("Philips .data file","*.data")])
outDir = Path(dlName).parent.absolute()
rlsName = filedialog.askopenfilename(title='Select .raw file',filetypes=[("Philips .raw file","*.raw")],initialdir=outDir)
path = os.path.normpath(rlsName)
fname = path.split(os.sep)

# Run converter
inputData = p2m.Ph2Mrd(dlName, rlsName)
inputData.trajtype = 2
inputData.delay = +2.5
mrdName, rls, dl = inputData.convert(outDir)
dset = mrd.Dataset(mrdName, "dataset", create_if_needed=False)

# Modify dset
header = mrd.xsd.CreateFromDocument(dset.read_xml_header())
studyInfo = mrd.xsd.studyInformationType()
subjectInfo = mrd.xsd.subjectInformationType()
enc = header.encoding[0]
pars = header.sequenceParameters
exp = header.experimentalConditions
sysInfo = header.acquisitionSystemInformation
trajDescr = mrd.xsd.trajectoryDescriptionType()
userParams = mrd.xsd.userParametersType()

sysInfo.institutionName = 'CCHMC' # hard code
exp.H1resonanceFrequency_Hz = 127753955 # hard code for 3T
sd = rls.header['sin']['start_scan_date_time'][0][0]
studyInfo.studyDate = datetime.datetime.strptime(sd,'%d-%b-%Y').strftime('%Y-%m-%d')
subjectInfo.patientID = fname[2]

if enc.trajectory == mrd.xsd.trajectoryType.RADIAL and float(rls.header['sin']['acq_gamma'][0][0]) < 42000.0: # set true flip angle for Xe acquisition
    pars.flipAngle_deg.insert(0,0.5)
    pars.flipAngle_deg.insert(1,20.0)


if 'Polarean_calibration'.lower() in rls.header['sin']['scan_name'][0][0].lower():
    pars.TR.insert(0,float(rls.header['sin']['repetition_times'][0][0]))
    pars.TR.insert(1,float(rls.header['sin']['repetition_times'][0][0]))
    pars.flipAngle_deg.insert(0,float(rls.header['sin']['flip_angles'][0][0]))
    pars.flipAngle_deg.insert(1,float(rls.header['sin']['flip_angles'][0][0]))
    

if 'Polarean_1ptDixon'.lower() in rls.header['sin']['scan_name'][0][0].lower():
    pars.TR.insert(0,float(rls.header['sin']['repetition_times'][0][0]) * 2)
    pars.TR.insert(1,float(rls.header['sin']['repetition_times'][0][0]) * 2)

if float(rls.header['sin']['acq_gamma'][0][0]) < 42000.0: 
    dissFreq = mrd.xsd.userParameterDoubleType('dissFreq_ppm')
    dissFreq.value = float(208)   
    userParams.userParameterDouble.insert(0,dissFreq)
else:
    sysInfo.systemFieldStrength_T = 3.0

orientation = mrd.xsd.userParameterStringType('orientation')
orientation.value = 'Coronal'   
userParams.userParameterString.insert(0,orientation)

centFreq = mrd.xsd.userParameterLongType('xe_center_frequency')
centFreq.value = float(rls.header['sin']['acq_gamma'][0][0]) / 42577.4688 * float(exp.H1resonanceFrequency_Hz)
userParams.userParameterLong.insert(0, centFreq)
offFreq = mrd.xsd.userParameterLongType('xe_dissolved_offset_frequency')
offFreq.value = float(rls.header['sin']['acq_gamma'][0][0]) / 42577.4688 * float(exp.H1resonanceFrequency_Hz) * 208 / 1000000
userParams.userParameterLong.insert(0, offFreq)


dwell = mrd.xsd.userParameterDoubleType('dwell')
dwell.value = float(rls.header['sin']['sample_time_interval'][0][0])
trajDescr.userParameterDouble.insert(0, dwell)

ramp_time = mrd.xsd.userParameterLongType('ramp_time')
try:
    ramp_time.value = round(float(rls.header['sin']['non_cart_fid_slope'][0][0])*dwell.value)
    trajDescr.userParameterLong.insert(0, ramp_time)
except:
    pass

header.userParameters = userParams
header.sequenceParameters = pars
header.acquisitionSystemInformation = sysInfo
header.experimentalConditions = exp
header.encoding[0].trajectoryDescription = trajDescr
header.studyInformation = studyInfo
header.subjectInformation = subjectInfo

# account for switching of data labels
if 'Polarean_calibration'.lower() in rls.header['sin']['scan_name'][0][0].lower():
    header.encoding[0].encodingLimits.contrast.minimum = 1
    header.encoding[0].encodingLimits.contrast.maximum = 2
    header.encoding[0].encodingLimits.contrast.center = 1
if 'Polarean_1ptDixon'.lower() in rls.header['sin']['scan_name'][0][0].lower():
    header.encoding[0].encodingLimits.contrast.minimum = 1
    header.encoding[0].encodingLimits.contrast.maximum = 2
    header.encoding[0].encodingLimits.contrast.center = 1
    header.encoding[0].encodingLimits.repetition.minimum = 0
    header.encoding[0].encodingLimits.repetition.maximum = 0
    header.encoding[0].encodingLimits.repetition.center = 0
    header.encoding[0].encodingLimits.set.minimum = 1
    header.encoding[0].encodingLimits.set.maximum = 2
    header.encoding[0].encodingLimits.set.center = 1

# account for 2nd echo in 1pt dixon
if 'Polarean_1ptDixon'.lower() in rls.header['sin']['scan_name'][0][0].lower():
    header.encoding.append(copy.deepcopy(header.encoding[0]))
    header.encoding[1].encodingLimits.kspace_encoding_step_0.minimum = -(header.encoding[1].encodingLimits.kspace_encoding_step_0.maximum+1)    
    header.encoding[0].encodedSpace.matrixSize.x = int(header.encoding[0].encodedSpace.matrixSize.x / 2)

dset.write_xml_header(mrd.xsd.ToXML(header))

for acqnum in range(dset.number_of_acquisitions()):
    acq_temp = dset.read_acquisition(acqnum)
    if 'Polarean_calibration'.lower() in rls.header['sin']['scan_name'][0][0].lower():
        # First 500 are dissolved; then gas
        if acqnum < 500:
            acq_temp.idx.contrast = 2
        elif acqnum >= 500:
            acq_temp.idx.contrast = 1
    if 'Polarean_1ptDixon'.lower() in rls.header['sin']['scan_name'][0][0].lower():
        # philips echoes = mrd contrasts = xemrd sets
        if acq_temp.idx.contrast == 0:
            acq_temp.idx.set = 1
        elif acq_temp.idx.contrast == 1:
            acq_temp.idx.set = 2
        # philips dynamics = mrd repititions = xemrd contrasts
            # 1 = gas/1; dynamic 2 = diss/2
        if acq_temp.idx.repetition == 0:
            acq_temp.idx.contrast = 1
        elif acq_temp.idx.repetition == 1:
            acq_temp.idx.contrast = 2
            acq_temp.idx.repetition = 0
    if 'Polarean_GX_UTE'.lower() in rls.header['sin']['scan_name'][0][0].lower():
        acq_temp.idx.contrast = 0
    # Replace old acq header
    dset.write_acquisition(acq_temp,acqnum)

# Save dset
dset.close()

# Rename
if 'Polarean_calibration'.lower() in rls.header['sin']['scan_name'][0][0].lower():
    try:
        os.remove(os.path.join(mrdName.parent, fname[2]+'_calibration.h5'))
    except:
        pass
    os.rename(mrdName, os.path.join(mrdName.parent, fname[2]+'_calibration.h5'))
if 'Polarean_1ptDixon'.lower() in rls.header['sin']['scan_name'][0][0].lower():
    try:
        os.remove(os.path.join(mrdName.parent, fname[2]+'_dixon.h5'))
    except:
        pass
    os.rename(mrdName, os.path.join(mrdName.parent, fname[2]+'_dixon.h5'))
if 'Polarean_GX_UTE'.lower() in rls.header['sin']['scan_name'][0][0].lower():
    try:
        os.remove(os.path.join(mrdName.parent, fname[2]+'_proton.h5'))
    except:
        pass
    os.rename(mrdName, os.path.join(mrdName.parent, fname[2]+'_proton.h5'))