from pathlib import Path
import sys
import os
import math
import ismrmrd as mrd
from xsdata.models.datatype import XmlDate, XmlTime

# OS Check
if sys.platform == "win32":
    pass
elif sys.platform == "darwin":
    raise RuntimeError('ReadPhilips not compiled for Mac')
elif sys.platform == "linux":
    raise RuntimeError('ReadPhilips not compiled for Linux')
else:
    raise RuntimeError('ReadPhilips not compiled for for OS')

# Python Version Check
if sys.version_info.major != 3:
    raise RuntimeError('Requiers python 3')
if sys.version_info.minor == 10:
    from rp.rp310win import *
elif sys.version_info.minor == 9:
    from rp.rp309win import *
elif sys.version_info.minor == 8:
    from rp.rp308win import *
elif sys.version_info.minor == 7:
    from rp.rp307win import *
elif sys.version_info.minor == 6:
    from rp.rp306win import *
else:
    raise RuntimeError('ReadPhilips not compiled for this python version')

class Ph2Mrd():
    def __init__(self, dlName, rlsName):               
        # check that the path actually exists
        if os.path.exists(dlName):
            self.dlName = dlName
        else:
            self.dlName = Path('')
        if os.path.exists(rlsName):
            self.rlsName = rlsName
        else:
            self.rlsName = Path('')
        
        # parameters for conversion
        self.trajtype = 0 # Use for user added orderings of trajectories
        self.delay = math.nan # Use for manual gr delays

    def convert(self, outDir):

        # Params
        dlPresent = False
        rlsPresent = False
        mrdName = 'philips_mrd.h5'

        # Process inputs
        dlFileName = Path(self.dlName)
        rlsFileName = Path(self.rlsName)
        outDir = Path(outDir)
        if dlFileName.exists():
            dlExt = dlFileName.suffix
            if str(dlExt).lower() in ['.list', '.data']:
                mrdName = dlFileName.stem
                mrdName = mrdName + '.h5'
                dlPresent = True
        if rlsFileName.exists():
            rlsExt = rlsFileName.suffix
            if str(rlsExt).lower() in ['.lab', '.raw', '.sin']:
                mrdName = rlsFileName.stem
                mrdName = mrdName + '.h5'
                rlsPresent = True
        if not Path(outDir).exists():
            outDir = Path().absolute()

        # Determine File Path
        mrdFileName = Path(outDir, mrdName)

        # Check that at least one file
        if not dlPresent and not rlsPresent:
            print('.data/.list and .raw/.lab/.sin not found.')
            raise FileNotFoundError('.data/.list and .raw/.lab/.sin not found.')
        elif dlPresent and not rlsPresent:
            print('.raw/.lab/.sin not found. Some header info will be missing.')
        elif not dlPresent and rlsPresent:
            print('.data/.list not found. Data may be uncorrected.')

        # Read in Philips data
        if dlPresent:
            dlPhData = PhilipsData(dlFileName)
            dlPhData.compute()
        if rlsPresent:
            rlsPhData = PhilipsData(rlsFileName)
            rlsPhData.trajtype = self.trajtype
            rlsPhData.delay = self.delay
            if dlPresent:
                rlsPhData.readParamOnly = True # use corrected data
            rlsPhData.compute()

        #FROM PhilipsData: outshape_string = np.array(['ch', 'mix', 'dyn', 'card', 'ex1', 'ex2',
        #                                               'echo', 'meas', 'loc', 'kz', 'ky', 'samp'])
        # Extract from philips data
        numChan = dlPhData.data.shape[0]
        numMix = dlPhData.data.shape[1]
        numDyn = dlPhData.data.shape[2]
        numCard = dlPhData.data.shape[3]
        numExtr1 = dlPhData.data.shape[4]
        numExtr2 = dlPhData.data.shape[5]
        numEcho = dlPhData.data.shape[6]
        numAver = dlPhData.data.shape[7]
        numLoca = dlPhData.data.shape[8]
        numKz = dlPhData.data.shape[9]
        numKy = dlPhData.data.shape[10]
        numKx = dlPhData.data.shape[11]

        dims = int(rlsPhData.header['sin']['encoding_dimensions'][0][0])
        try:
            traj_type = int(rlsPhData.header['sin']['k_space_traj_type'][0][0])
        except:
            traj_type = 0
        if traj_type == 1:  # radial
            crds = rlsPhData.radparams['COORDS']
            min_samp = int(rlsPhData.header['sin']['non_cart_min_encoding_nrs'][0][0])
            cent_samp = 0
            max_samp = int(rlsPhData.header['sin']['non_cart_max_encoding_nrs'][0][0])
            min_ky = int(rlsPhData.header['sin']['non_cart_min_encoding_nrs'][0][1])
            cent_ky = 0
            max_ky = int(rlsPhData.header['sin']['non_cart_max_encoding_nrs'][0][1])
            min_kz = int(rlsPhData.header['sin']['non_cart_min_encoding_nrs'][0][2])
            cent_kz = 0
            max_kz = int(rlsPhData.header['sin']['non_cart_max_encoding_nrs'][0][2])
        elif traj_type == 2:  # spiral
            crds = rlsPhData.spparams['COORDS_EXPANDED']
            min_samp = int(rlsPhData.header['sin']['non_cart_min_encoding_nrs'][0][0])
            cent_samp = 0
            max_samp = int(rlsPhData.header['sin']['non_cart_max_encoding_nrs'][0][0])
            min_ky = int(rlsPhData.header['sin']['non_cart_min_encoding_nrs'][0][1])
            cent_ky = 0
            max_ky = int(rlsPhData.header['sin']['non_cart_max_encoding_nrs'][0][1])
            min_kz = int(rlsPhData.header['sin']['non_cart_min_encoding_nrs'][0][2])
            cent_kz = 0
            max_kz = int(rlsPhData.header['sin']['non_cart_max_encoding_nrs'][0][2])
        else:
            min_samp = int(rlsPhData.header['sin']['min_encoding_numbers'][0][0])
            cent_samp = 0
            max_samp = int(rlsPhData.header['sin']['max_encoding_numbers'][0][0])
            min_ky = int(rlsPhData.header['sin']['min_encoding_numbers'][0][1])
            cent_ky = 0
            max_ky = int(rlsPhData.header['sin']['max_encoding_numbers'][0][1])
            min_kz = int(rlsPhData.header['sin']['min_encoding_numbers'][0][2])
            cent_kz = 0
            max_kz = int(rlsPhData.header['sin']['max_encoding_numbers'][0][2])

        # Open the dataset
        if mrdFileName.exists():
            mrdFileName.unlink() # delete if already exists
        dset = mrd.Dataset(mrdFileName, "dataset", create_if_needed=True)
        
        # Create the XML header and write it to the file
        header = mrd.xsd.ismrmrdHeader()
        
        # Experimental Conditions
        exp = mrd.xsd.experimentalConditionsType()
        exp.H1resonanceFrequency_Hz = 127728000
        header.experimentalConditions = exp
        
        # Acquisition System Information
        sys = mrd.xsd.acquisitionSystemInformationType()
        sys.systemVendor = 'Philips'
        sys.receiverChannels = numChan
        if float(rlsPhData.header['sin']['acq_gamma'][0][0]) < 42000.0: #MN only offered on 3T
            sys.systemFieldStrength_T = 3.0
        header.acquisitionSystemInformation = sys

        # Measurement Information
        meas_info = mrd.xsd.measurementInformationType()
        scan_date = mrdName[:8]
        scan_time = mrdName[9:15]
        meas_info.frameOfReferenceUID = scan_date
        meas_info.protocolName = rlsPhData.header['sin']['scan_name'][0][0]
        meas_info.seriesDate = XmlDate(int(scan_date[:4]), int(scan_date[4:6]), int(scan_date[6:]))
        meas_info.seriesTime = XmlTime(int(scan_time[:2]), int(scan_time[2:4]), int(scan_time[4:6]))
        header.measurementInformation = meas_info

        # Sequence Parameters
        pars = mrd.xsd.sequenceParametersType()
        pars.TE.insert(0,float(rlsPhData.header['sin']['echo_times'][0][0]))
        pars.TI.insert(0,float(rlsPhData.header['sin']['inversion_delays'][0][0]))
        pars.TR.insert(0,float(rlsPhData.header['sin']['repetition_times'][0][0]))
        if numEcho > 1:
            pars.echo_spacing.insert(0,float(rlsPhData.header['sin']['echo_times'][1][0])-float(rlsPhData.header['sin']['echo_times'][0][0]))
        pars.flipAngle_deg.insert(0,float(rlsPhData.header['sin']['flip_angles'][0][0]))
        header.sequenceParameters = pars
        
        # Encoding
        encoding = mrd.xsd.encodingType()
        if traj_type == 0:
            encoding.trajectory = mrd.xsd.trajectoryType.CARTESIAN     
        if traj_type == 1:
            encoding.trajectory = mrd.xsd.trajectoryType.RADIAL 
        if traj_type == 2:
            encoding.trajectory = mrd.xsd.trajectoryType.SPIRAL
        
        # encoded and recon spaces; assuming no change in FOV between encoded and reconned
        rfov = mrd.xsd.fieldOfViewMm()
        rfov.x = float(rlsPhData.header['sin']['recon_resolutions'][0][0]) * float(rlsPhData.header['sin']['voxel_sizes'][0][0])
        rfov.y = float(rlsPhData.header['sin']['recon_resolutions'][0][1]) * float(rlsPhData.header['sin']['voxel_sizes'][0][1])
        rfov.z = float(rlsPhData.header['sin']['recon_resolutions'][0][2]) * float(rlsPhData.header['sin']['voxel_sizes'][0][2])
        efov = mrd.xsd.fieldOfViewMm()
        efov.x = float(rlsPhData.header['sin']['oversample_factors'][0][0]) * float(rfov.x)
        efov.y = float(rlsPhData.header['sin']['oversample_factors'][0][1]) * float(rfov.y)
        efov.z = float(rlsPhData.header['sin']['oversample_factors'][0][2]) * float(rfov.z)

        ematrix = mrd.xsd.matrixSizeType()
        ematrix.x = numKx
        ematrix.y = numKy
        ematrix.z = numKz
        rmatrix = mrd.xsd.matrixSizeType()
        rmatrix.x = int(rlsPhData.header['sin']['recon_resolutions'][0][0])
        rmatrix.y = int(rlsPhData.header['sin']['recon_resolutions'][0][1])
        rmatrix.z = int(rlsPhData.header['sin']['recon_resolutions'][0][2])
        
        espace = mrd.xsd.encodingSpaceType()
        espace.matrixSize = ematrix
        espace.fieldOfView_mm = efov
        rspace = mrd.xsd.encodingSpaceType()
        rspace.matrixSize = rmatrix
        rspace.fieldOfView_mm = rfov
        
        # Set encoded and recon spaces
        encoding.encodedSpace = espace
        encoding.reconSpace = rspace
        
        # Encoding limits
        limits = mrd.xsd.encodingLimitsType()
        
        limitsA = mrd.xsd.limitType()
        limitsA.minimum = 0
        limitsA.center = round(numAver / 2)
        limitsA.maximum = numAver - 1
        limits.average = limitsA

        limitsC = mrd.xsd.limitType()
        limitsC.minimum = 0
        limitsC.center = round(numEcho / 2)
        limitsC.maximum = numEcho - 1
        limits.contrast = limitsC

        limits1 = mrd.xsd.limitType()
        limits1.minimum = min_ky
        limits1.center = cent_ky
        limits1.maximum = max_ky
        limits.kspace_encoding_step_1 = limits1
        
        limits2 = mrd.xsd.limitType()
        limits2.minimum = min_kz
        limits2.center = cent_kz
        limits2.maximum = max_kz
        limits.kspace_encoding_step_2 = limits2

        limitsP = mrd.xsd.limitType()
        limitsP.minimum = 0
        limitsP.center = round(numCard / 2)
        limitsP.maximum = numCard - 1
        limits.phase = limitsP

        limitsR = mrd.xsd.limitType()
        limitsR.minimum = 0
        limitsR.center = round(numDyn / 2)
        limitsR.maximum = numDyn - 1
        limits.repetition = limitsR

        limitsEx1 = mrd.xsd.limitType()
        limitsEx1.minimum = 0
        limitsEx1.center = round(numExtr1 / 2)
        limitsEx1.maximum = numExtr1 - 1
        limits.segment = limitsEx1

        limitsSet = mrd.xsd.limitType()
        limitsSet.minimum = 0
        limitsSet.center = round(numMix / 2)
        limitsSet.maximum = numMix - 1
        limits.set = limitsSet

        limitsSl = mrd.xsd.limitType()
        limitsSl.minimum = 0
        limitsSl.center = round(numLoca / 2)
        limitsSl.maximum = numLoca - 1
        limits.slice = limitsSl
        
        limits0 = mrd.xsd.limitType()
        limits0.minimum = min_samp
        limits0.center = cent_samp
        limits0.maximum = max_samp
        limits.kspace_encoding_step_0 = limits0
        
        encoding.encodingLimits = limits
        
        # append encoding to header
        header.encoding.append(encoding)

        # write header
        dset.write_xml_header(mrd.xsd.ToXML(header))

        # add data
        acq = mrd.Acquisition()
        acq_head = mrd.AcquisitionHeader()
        acq_head.number_of_samples = numKx
        acq_head.active_channels = numChan
        acq_head.trajectory_dimensions = dims
        acq_head.sample_time_us = float(rlsPhData.header['sin']['sample_time_interval'][0][0])
        acq.resize(numKx, numChan)
        acq.version = 1
        acq.available_channels = numChan
        acq.center_sample = cent_samp
        acq.read_dir[0] = 1.0
        acq.phase_dir[1] = 1.0
        acq.slice_dir[2] = 1.0
        acq.setHead(acq_head)

        for a in range(dlPhData.header['list']['typ'].size):
            # Reset
            acq.clear_all_flags()

            # Index
            acq.scan_counter = a
            acq.idx.average = int(dlPhData.header['list']['aver'][a])
            acq.idx.contrast = int(dlPhData.header['list']['echo'][a])
            acq.idx.kspace_encode_step_1 = int(dlPhData.header['list']['ky'][a]) - min_ky
            acq.idx.kspace_encode_step_2 = int(dlPhData.header['list']['kz'][a]) - min_kz
            acq.idx.phase = int(dlPhData.header['list']['card'][a])
            acq.idx.repetition = int(dlPhData.header['list']['dyn'][a])
            acq.idx.segment = int(dlPhData.header['list']['extr1'][a])
            acq.idx.set = int(dlPhData.header['list']['mix'][a])
            acq.idx.slice = int(dlPhData.header['list']['loca'][a])

            if dlPhData.header['list']['typ'][a] != 'STD':
                continue #only worrying about std data for now
            if int(dlPhData.header['list']['extr2'][a]) > 0:
                continue #dont have dimension for it            
            if int(dlPhData.header['list']['chan'][a]) > 0:
                continue # already added by chan 0

            # Data
            #FROM PhilipsData: outshape_string = np.array(['ch', 'mix', 'dyn', 'card', 'ex1', 'ex2',
            #                                               'echo', 'meas', 'loc', 'kz', 'ky', 'samp'])
            dat = dlPhData.data[:,#'ch'
            acq.idx.set,#'mix'
            acq.idx.repetition,#'dyn',
            acq.idx.phase,#'card',
            acq.idx.segment,#'ex1',
            int(dlPhData.header['list']['extr2'][a]),#'ex2',
            acq.idx.contrast,#'echo',
            acq.idx.average,#'meas',
            acq.idx.slice,#'loc',
            acq.idx.kspace_encode_step_2,#'kz',
            acq.idx.kspace_encode_step_1,#'ky',
            :]#'samp'
            acq.data[:] = dat
            try:
                traj = crds[acq.idx.kspace_encode_step_2, acq.idx.kspace_encode_step_1, :, :]
                acq.traj[:] = traj
            except:
                pass

            # Flags
            if acq.idx.repetition == 0:
                acq.setFlag(mrd.ACQ_FIRST_IN_REPETITION)
            elif acq.idx.repetition == numDyn - 1:
                acq.setFlag(mrd.ACQ_LAST_IN_REPETITION)
            if acq.idx.kspace_encode_step_1 == 0:
                acq.setFlag(mrd.ACQ_FIRST_IN_ENCODE_STEP1)
            elif acq.idx.kspace_encode_step_1 == numKy - 1:
                acq.setFlag(mrd.ACQ_LAST_IN_ENCODE_STEP1)
            if acq.idx.kspace_encode_step_2 == 0:
                acq.setFlag(mrd.ACQ_FIRST_IN_ENCODE_STEP2)
            elif acq.idx.kspace_encode_step_2 == numKz - 1:
                acq.setFlag(mrd.ACQ_LAST_IN_ENCODE_STEP2)

            # Add
            dset.append_acquisition(acq)    
        
        dset.close()

        return mrdFileName, rlsPhData, dlPhData
