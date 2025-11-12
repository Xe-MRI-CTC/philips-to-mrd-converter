from pathlib import Path
import sys
import os
import math
import importlib
import ismrmrd as mrd
from xsdata.models.datatype import XmlDate, XmlTime

# Import read philips
if sys.version_info.major != 3:
    raise RuntimeError('Requires python 3')

major_version = sys.version_info.major
minor_version = sys.version_info.minor
rp_name = f"rp.rp{major_version}{minor_version}"
try:
    rp = importlib.import_module(rp_name)
except ModuleNotFoundError:
    raise RuntimeError(
        f'ReadPhilips not compiled for Python {sys.version_info.major}.{sys.version_info.minor}')


class Ph2Mrd():
    def __init__(self, dlName=None, rlsName=None):
        # check that the path actually exists
        if dlName != None and os.path.exists(dlName):
            self.dlName = dlName
        else:
            self.dlName = Path('')
        if rlsName != None and os.path.exists(rlsName):
            self.rlsName = rlsName
        else:
            self.rlsName = Path('')

        # parameters for conversion
        self.trajorder = 0  # Use for user added orderings of trajectories
        self.delay = math.nan  # Use for manual gr delays

    def convert(self, outDir):

        # Default Params
        mrdName = 'philips_mrd.h5'

        # Process inputs
        outDir, dlPresent, rlsPresent, mrdName, dlFileName, rlsFileName = self.get_inputs(
            outDir)

        # Check that at least one raw data type provided
        if not dlPresent and not rlsPresent:
            print('.data/.list and .raw/.lab/.sin not found.')
            raise FileNotFoundError(
                '.data/.list and/or .raw/.lab/.sin not found.')

        # Determine File Path
        mrdFileName = Path(outDir, mrdName)

        # Read in Philips data
        dlPhData, rlsPhData = self.read_data(
            dlPresent, rlsPresent, dlFileName, rlsFileName)

        # Extract from philips data size
        data_size, traj_type = self.get_data_sizes(
            dlPresent, dlPhData, rlsPresent, rlsPhData)

        # Extract Non cartesian sizes/trajectory
        crds, crds_flyback = self.get_non_cart_coords(
            rlsPhData, traj_type, data_size)

        # Open the dataset
        if mrdFileName.exists():
            mrdFileName.unlink()  # delete if already exists
        dset = mrd.Dataset(mrdFileName, "dataset", create_if_needed=True)

        # Create the XML header and write it to the file
        header = mrd.xsd.ismrmrdHeader()

        # Experimental Conditions
        self.setExperiment(header)

        # Acquisition System Information
        self.setSystem(rlsPhData, data_size, header)

        # Measurement Information
        self.setMeasurementInfo(mrdName, rlsPhData, header)

        # Sequence Parameters
        self.setSequenceParameters(rlsPhData, data_size, header)

        # Encoding
        self.setEncoding(rlsPhData, data_size, traj_type, header)

        # write header
        dset.write_xml_header(mrd.xsd.ToXML(header))

        # add data
        self.addAcquisitions(dlPhData, rlsPhData, data_size,
                             crds, crds_flyback, dset)

        return mrdFileName, rlsPhData, dlPhData

    def addAcquisitions(self, dlPhData, rlsPhData, data_size, crds, crds_flyback, dset):
        acq = mrd.Acquisition()
        acq_head = mrd.AcquisitionHeader()
        acq_head.number_of_samples = data_size.numKx
        acq_head.active_channels = data_size.numChan
        acq_head.trajectory_dimensions = data_size.dims
        try:
            acq_head.sample_time_us = float(
                rlsPhData.header['sin']['sample_time_interval'][0][0])
        except:
            # don't know so make t-axis = points
            acq_head.sample_time_us = float(1.0)
            print('Warning: dwell time not found; setting to 1us.')
        acq.version = 1
        acq.available_channels = data_size.numChan
        acq.center_sample = data_size.kxCent
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
            acq.idx.kspace_encode_step_1 = int(
                dlPhData.header['list']['ky'][a]) - data_size.kyMin
            acq.idx.kspace_encode_step_2 = int(
                dlPhData.header['list']['kz'][a]) - data_size.kzMin
            acq.idx.phase = int(dlPhData.header['list']['card'][a])
            acq.idx.repetition = int(dlPhData.header['list']['dyn'][a])
            acq.idx.segment = int(dlPhData.header['list']['extr1'][a])
            acq.idx.set = int(dlPhData.header['list']['mix'][a])
            acq.idx.slice = int(dlPhData.header['list']['loca'][a])

            if dlPhData.header['list']['typ'][a] != 'STD':
                continue  # only worrying about std data for now
            if int(dlPhData.header['list']['extr2'][a]) > 0:
                continue  # dont have dimension for it
            if int(dlPhData.header['list']['chan'][a]) > 0:
                continue  # already added by chan 0

            # Set up data
            if acq.idx.contrast == 0:
                nKx = data_size.numKx
            else:
                nKx = data_size.numKxEcho
            acq.resize(nKx, data_size.numChan, data_size.dims)

            # Data
            # FROM PhilipsData: outshape_string = np.array(['ch', 'mix', 'dyn', 'card', 'ex1', 'ex2',
            #                                               'echo', 'meas', 'loc', 'kz', 'ky', 'samp'])
            dat = dlPhData.data[:,  # 'ch'
                                acq.idx.set,  # 'mix'
                                acq.idx.repetition,  # 'dyn',
                                acq.idx.phase,  # 'card',
                                acq.idx.segment,  # 'ex1',
                                int(dlPhData.header['list']
                                    ['extr2'][a]),  # 'ex2',
                                acq.idx.contrast,  # 'echo',
                                acq.idx.average,  # 'meas',
                                acq.idx.slice,  # 'loc',
                                acq.idx.kspace_encode_step_2,  # 'kz',
                                acq.idx.kspace_encode_step_1,  # 'ky',
                                :]  # 'samp'
            acq.data[:] = dat[:, :nKx]
            try:
                if acq.idx.contrast == 0:
                    traj = crds[acq.idx.kspace_encode_step_2,
                                acq.idx.kspace_encode_step_1, :, :]
                else:
                    traj = crds_flyback[acq.idx.kspace_encode_step_2,
                                        acq.idx.kspace_encode_step_1, :, :]
                acq.traj[:] = traj
            except:
                pass

            # Flags
            if acq.idx.repetition == 0:
                acq.setFlag(mrd.ACQ_FIRST_IN_REPETITION)
            elif acq.idx.repetition == data_size.numDyn - 1:
                acq.setFlag(mrd.ACQ_LAST_IN_REPETITION)
            if acq.idx.kspace_encode_step_1 == 0:
                acq.setFlag(mrd.ACQ_FIRST_IN_ENCODE_STEP1)
            elif acq.idx.kspace_encode_step_1 == data_size.numKy - 1:
                acq.setFlag(mrd.ACQ_LAST_IN_ENCODE_STEP1)
            if acq.idx.kspace_encode_step_2 == 0:
                acq.setFlag(mrd.ACQ_FIRST_IN_ENCODE_STEP2)
            elif acq.idx.kspace_encode_step_2 == data_size.numKz - 1:
                acq.setFlag(mrd.ACQ_LAST_IN_ENCODE_STEP2)

            # Add
            dset.append_acquisition(acq)

        dset.close()

    def setEncoding(self, rlsPhData, data_size, traj_type, header):
        encoding = mrd.xsd.encodingType()
        if traj_type == 0:
            encoding.trajectory = mrd.xsd.trajectoryType.CARTESIAN
        if traj_type == 1:
            encoding.trajectory = mrd.xsd.trajectoryType.RADIAL
        if traj_type == 2:
            encoding.trajectory = mrd.xsd.trajectoryType.SPIRAL

        # encoded and recon spaces; assuming no change in FOV between encoded and reconned
        self.setEncodingSpace(rlsPhData, data_size, encoding)

        # Encoding limits
        self.setEncodingLimits(data_size, encoding)

        # append encoding to header
        header.encoding.append(encoding)

    def setEncodingLimits(self, data_size, encoding):
        limits = mrd.xsd.encodingLimitsType()

        limitsA = mrd.xsd.limitType()
        limitsA.minimum = 0
        limitsA.center = round(data_size.numAver / 2)
        limitsA.maximum = data_size.numAver - 1
        limits.average = limitsA

        limitsC = mrd.xsd.limitType()
        limitsC.minimum = 0
        limitsC.center = round(data_size.numEcho / 2)
        limitsC.maximum = data_size.numEcho - 1
        limits.contrast = limitsC

        limits1 = mrd.xsd.limitType()
        limits1.minimum = data_size.kyMin
        limits1.center = data_size.kyCent
        limits1.maximum = data_size.kyMax
        limits.kspace_encoding_step_1 = limits1

        limits2 = mrd.xsd.limitType()
        limits2.minimum = data_size.kzMin
        limits2.center = data_size.kzCent
        limits2.maximum = data_size.kzMax
        limits.kspace_encoding_step_2 = limits2

        limitsP = mrd.xsd.limitType()
        limitsP.minimum = 0
        limitsP.center = round(data_size.numCard / 2)
        limitsP.maximum = data_size.numCard - 1
        limits.phase = limitsP

        limitsR = mrd.xsd.limitType()
        limitsR.minimum = 0
        limitsR.center = round(data_size.numDyn / 2)
        limitsR.maximum = data_size.numDyn - 1
        limits.repetition = limitsR

        limitsEx1 = mrd.xsd.limitType()
        limitsEx1.minimum = 0
        limitsEx1.center = round(data_size.numExtr1 / 2)
        limitsEx1.maximum = data_size.numExtr1 - 1
        limits.segment = limitsEx1

        limitsSet = mrd.xsd.limitType()
        limitsSet.minimum = 0
        limitsSet.center = round(data_size.numMix / 2)
        limitsSet.maximum = data_size.numMix - 1
        limits.set = limitsSet

        limitsSl = mrd.xsd.limitType()
        limitsSl.minimum = 0
        limitsSl.center = round(data_size.numLoca / 2)
        limitsSl.maximum = data_size.numLoca - 1
        limits.slice = limitsSl

        limits0 = mrd.xsd.limitType()
        limits0.minimum = data_size.kxMin
        limits0.center = data_size.kxCent
        limits0.maximum = data_size.kxMax
        limits.kspace_encoding_step_0 = limits0

        encoding.encodingLimits = limits

    def setEncodingSpace(self, rlsPhData, data_size, encoding):
        rfov = mrd.xsd.fieldOfViewMm()
        rfov.x = float(rlsPhData.header['sin']['recon_resolutions'][0]
                       [0]) * float(rlsPhData.header['sin']['voxel_sizes'][0][0])
        rfov.y = float(rlsPhData.header['sin']['recon_resolutions'][0]
                       [1]) * float(rlsPhData.header['sin']['voxel_sizes'][0][1])
        rfov.z = float(rlsPhData.header['sin']['recon_resolutions'][0]
                       [2]) * float(rlsPhData.header['sin']['voxel_sizes'][0][2])
        efov = mrd.xsd.fieldOfViewMm()
        efov.x = float(
            rlsPhData.header['sin']['oversample_factors'][0][0]) * float(rfov.x)
        efov.y = float(
            rlsPhData.header['sin']['oversample_factors'][0][1]) * float(rfov.y)
        efov.z = float(
            rlsPhData.header['sin']['oversample_factors'][0][2]) * float(rfov.z)

        ematrix = mrd.xsd.matrixSizeType()
        ematrix.x = data_size.numKx
        ematrix.y = data_size.numKy
        ematrix.z = data_size.numKz
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

    def setSequenceParameters(self, rlsPhData, data_size, header):
        pars = mrd.xsd.sequenceParametersType()
        pars.TE.insert(0, float(rlsPhData.header['sin']['echo_times'][0][0]))
        pars.TI.insert(
            0, float(rlsPhData.header['sin']['inversion_delays'][0][0]))
        pars.TR.insert(
            0, float(rlsPhData.header['sin']['repetition_times'][0][0]))
        if data_size.numEcho > 1:
            pars.echo_spacing.insert(0, float(
                rlsPhData.header['sin']['echo_times'][0][1])-float(rlsPhData.header['sin']['echo_times'][0][0]))
        pars.flipAngle_deg.insert(
            0, float(rlsPhData.header['sin']['flip_angles'][0][0]))
        header.sequenceParameters = pars

    def setMeasurementInfo(self, mrdName, rlsPhData, header):
        meas_info = mrd.xsd.measurementInformationType()
        scan_date = mrdName[:8]
        scan_time = mrdName[9:15]
        meas_info.frameOfReferenceUID = scan_date
        meas_info.protocolName = rlsPhData.header['sin']['scan_name'][0][0]
        meas_info.seriesDate = XmlDate(
            int(scan_date[:4]), int(scan_date[4:6]), int(scan_date[6:]))
        meas_info.seriesTime = XmlTime(
            int(scan_time[:2]), int(scan_time[2:4]), int(scan_time[4:6]))
        header.measurementInformation = meas_info

    def setSystem(self, rlsPhData, data_size, header):
        sys = mrd.xsd.acquisitionSystemInformationType()
        sys.systemVendor = 'Philips'
        sys.receiverChannels = data_size.numChan
        # hard coded as MN only offered on 3T
        if float(rlsPhData.header['sin']['acq_gamma'][0][0]) < 42000.0:
            sys.systemFieldStrength_T = 3.0
        header.acquisitionSystemInformation = sys

    def setExperiment(self, header):
        exp = mrd.xsd.experimentalConditionsType()
        exp.H1resonanceFrequency_Hz = 127728000  # hard coded to 3T
        header.experimentalConditions = exp

    def get_non_cart_coords(self, rlsPhData, traj_type, data_size):
        crds = []
        crds_flyback = []
        if traj_type == 1:  # radial
            crds = rlsPhData.radparams['COORDS']
            if data_size.numEcho > 1:
                crds_flyback = rlsPhData.radparams['COORDS_FLYBACK']
        elif traj_type == 2:  # spiral
            crds = rlsPhData.spparams['COORDS_EXPANDED']
        return crds, crds_flyback

    def get_data_sizes(self, dlPresent, dlPhData, rlsPresent, rlsPhData):
        # Determine if Cart or non-Cart
        try:
            traj_type = int(rlsPhData.header['sin']['k_space_traj_type'][0][0])
        except:
            traj_type = 0
        # Create data size instance and update
        data_size = PhDataSize()
        raw_data_source = None
        # FROM PhilipsData: outshape_string = np.array(['ch', 'mix', 'dyn', 'card', 'ex1', 'ex2',
        #                                               'echo', 'meas', 'loc', 'kz', 'ky', 'samp'])
        if dlPresent:
            raw_data_source = dlPhData
            # TODO Needs imtrovment for accuracy
            data_size.dims = 3 if len(
                set(dlPhData.header['list']['kz'])) > 1 else 2
        elif rlsPresent:
            raw_data_source = rlsPhData
            data_size.dims = int(
                rlsPhData.header['sin']['encoding_dimensions'][0][0])
        data_size.numChan = raw_data_source.data.shape[0]
        data_size.numMix = raw_data_source.data.shape[1]
        data_size.numDyn = raw_data_source.data.shape[2]
        data_size.numCard = raw_data_source.data.shape[3]
        data_size.numExtr1 = raw_data_source.data.shape[4]
        data_size.numExtr2 = raw_data_source.data.shape[5]
        data_size.numEcho = raw_data_source.data.shape[6]
        data_size.numAver = raw_data_source.data.shape[7]
        data_size.numLoca = raw_data_source.data.shape[8]
        data_size.numKz = raw_data_source.data.shape[9]
        data_size.numKy = raw_data_source.data.shape[10]
        data_size.numKx = raw_data_source.data.shape[11]

        # for mEcho cases with flyback
        if dlPresent and len(set(dlPhData.header['list']['size'])) > 2:
            data_size.numKxEcho = data_size.numKx
            data_size.numKx = int(data_size.numKx / 2)
        # for mEcho cases with flyback
        elif rlsPresent and len(set(rlsPhData.header['sin']['data_burst_size'][0])) > 1:
            data_size.numKxEcho = data_size.numKx
            data_size.numKx = int(data_size.numKx / 2)

        enc_nr_name_min = 'min_encoding_numbers'
        enc_nr_name_max = 'max_encoding_numbers'
        if traj_type > 0:
            enc_nr_name_min = 'non_cart_min_encoding_nrs'
            enc_nr_name_max = 'non_cart_max_encoding_nrs'

        data_size.kxMin = int(rlsPhData.header['sin'][enc_nr_name_min][0][0])
        data_size.kxMax = int(rlsPhData.header['sin'][enc_nr_name_max][0][0])
        data_size.kxCent = 0
        data_size.kyMin = int(rlsPhData.header['sin'][enc_nr_name_min][0][1])
        data_size.kyMax = int(rlsPhData.header['sin'][enc_nr_name_max][0][1])
        data_size.kxCent = 0
        data_size.kzMin = int(rlsPhData.header['sin'][enc_nr_name_min][0][2])
        data_size.kzMax = int(rlsPhData.header['sin'][enc_nr_name_max][0][2])
        data_size.kzCent = 0

        return data_size, traj_type

    def read_data(self, dlPresent, rlsPresent, dlFileName, rlsFileName):
        dlPhData = None
        rlsPhData = None
        if dlPresent:
            dlPhData = rp.PhilipsData(dlFileName)
            dlPhData.compute()
        if rlsPresent:
            rlsPhData = rp.PhilipsData(rlsFileName)
            rlsPhData.trajtype = self.trajorder
            rlsPhData.delay = self.delay
            if dlPresent:
                rlsPhData.readParamOnly = True  # use corrected data
            rlsPhData.compute()
        return dlPhData, rlsPhData

    def get_inputs(self, outDir):
        dlPresent = False
        rlsPresent = False
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
        return outDir, dlPresent, rlsPresent, mrdName, dlFileName, rlsFileName


class PhDataSize():
    def __init__(self):
        self.numChan = 0
        self.numMix = 0
        self.numDyn = 0
        self.numCard = 0
        self.numExtr1 = 0
        self.numExtr2 = 0
        self.numEcho = 0
        self.numAver = 0
        self.numLoca = 0
        self.kzMin = 0
        self.kzMax = 0
        self.kzCent = 0
        self.numKz = 0
        self.kyMin = 0
        self.kyMax = 0
        self.kyCent = 0
        self.numKy = 0
        self.kxMin = 0
        self.kxMax = 0
        self.kxCent = 0
        self.numKx = 0
        self.numKxEcho = 0
        self.dims = 0
