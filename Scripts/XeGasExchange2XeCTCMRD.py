import sys  # noqa: E402
from pathlib import Path  # noqa: E402
p2mDir = Path(__file__).parent.parent.absolute()  # noqa: E402
sys.path.append(str(p2mDir))  # noqa: E402
import philips2mrd as p2m
import ismrmrd as mrd
from Scripts.XeGasExchange2XeCTCMRD_Config import Config, DataType
import argparse
from tkinter import filedialog
import tkinter as tk
import numpy as np
import copy
import os

# Constants
H1_GAMMA = 42577.4688


def Gx2XeCTCMRD(data_file=None, raw_file=None, traj_file=None):
    # Get paths
    root = tk.Tk()
    root.withdraw()
    if data_file == None:
        dlName = filedialog.askopenfilename(title='Select .data file', filetypes=[
            ("Philips .data file", "*.data")])
    else:
        dlName = data_file
    outDir = Path(dlName).parent.absolute()
    if raw_file == None:
        rlsName = filedialog.askopenfilename(title='Select .raw file', filetypes=[
            ("Philips .raw file", "*.raw")], initialdir=outDir)
    else:
        rlsName = raw_file
    path = os.path.normpath(rlsName)
    fname = path.split(os.sep)
    try:
        patientID = fname[-2]  # assume patient ID is name of folder
    except:
        print('Folder not found.')
        raise FileNotFoundError('Folder not found.')

    # Get config
    data_set_config = Config()

    # Run converter
    inputData = p2m.Ph2Mrd(dlName, rlsName)
    inputData.trajorder = data_set_config.trajorder
    inputData.delay = data_set_config.gr_delay
    mrdName, rls, dl = inputData.convert(outDir)
    dset = mrd.Dataset(mrdName, "dataset", create_if_needed=False)

    # Get Config for XeMRD scan
    header = mrd.xsd.CreateFromDocument(dset.read_xml_header())
    data_set_config.update(dl, rls, header)

    # Get path to sin file trajectories if necessary
    if data_set_config.ext_traj == True:
        if traj_file == None:
            ext_coords = filedialog.askopenfilename(title='Select trajectory .sin file', filetypes=[
                ("Trajectory .sin file", "*.sin")], initialdir=outDir)
        else:
            ext_coords = traj_file
        inputDataTraj = p2m.Ph2Mrd(dlName, ext_coords)
        inputDataTraj.trajorder = data_set_config.trajorder
        inputDataTraj.delay = data_set_config.gr_delay
        mrdNameCoords, rlsCoords, _ = inputDataTraj.convert(outDir)
        os.remove(mrdNameCoords)
        try:
            traj_type = int(rlsCoords.header['sin']['k_space_traj_type'][0][0])
        except:
            traj_type = 0
        if traj_type == 1:  # radial
            crds = rlsCoords.radparams['COORDS']
            if int(rlsCoords.header['sin']['nr_echoes'][0][0]) > 1:
                crds_flyback = rlsCoords.radparams['COORDS_FLYBACK']
        elif traj_type == 2:  # spiral
            crds = rlsCoords.spparams['COORDS_EXPANDED']

    # Get dset header
    studyInfo = header.studyInformation
    subjectInfo = mrd.xsd.subjectInformationType()
    enc = header.encoding[0]
    pars = header.sequenceParameters
    exp = header.experimentalConditions
    sysInfo = header.acquisitionSystemInformation
    trajDescr = mrd.xsd.trajectoryDescriptionType()
    userParams = mrd.xsd.userParametersType()

    # Modify basic info
    sysInfo.institutionName = data_set_config.institution
    sysInfo.systemFieldStrength_T = data_set_config.field_strength
    exp.H1resonanceFrequency_Hz = data_set_config.H1resonanceFrequency_Hz
    subjectInfo.patientID = patientID

    orientation = mrd.xsd.userParameterStringType('orientation')
    orientation.value = data_set_config.orientation
    userParams.userParameterString.insert(0, orientation)

    # Modify to xenon MRD header
    dwell = mrd.xsd.userParameterDoubleType('dwell')
    if data_set_config.ext_traj == True:
        dwell.value = float(
            rlsCoords.header['sin']['sample_time_interval'][0][0])
    else:
        dwell.value = float(rls.header['sin']['sample_time_interval'][0][0])
    trajDescr.userParameterDouble.insert(0, dwell)

    ramp_time = mrd.xsd.userParameterLongType('ramp_time')
    try:
        ramp_time.value = round(
            float(rls.header['sin']['non_cart_fid_slope'][0][0])*dwell.value)
        trajDescr.userParameterLong.insert(0, ramp_time)
    except:
        pass

    if data_set_config.data_type == DataType.DIXON:  # set true flip angle for Xe Dixon acqs
        pars.flipAngle_deg.insert(0, data_set_config.flip_angle_gas)
        pars.flipAngle_deg.insert(1, data_set_config.flip_angle_dis)

    if data_set_config.data_type != DataType.UTE:  # xenon has two trs for two frequencies
        # tr_factor accounts for different frequencies
        pars.TR.insert(0, float(
            rls.header['sin']['repetition_times'][0][0]) * data_set_config.tr_factor)
        pars.TR.insert(1, float(
            rls.header['sin']['repetition_times'][0][0]) * data_set_config.tr_factor)

    if data_set_config.data_type != DataType.UTE:  # xenon needs center freq and offset
        centFreq = mrd.xsd.userParameterLongType('xe_center_frequency')
        centFreq.value = float(
            rls.header['sin']['acq_gamma'][0][0]) / H1_GAMMA * float(exp.H1resonanceFrequency_Hz)
        userParams.userParameterLong.insert(0, centFreq)
        offFreq = mrd.xsd.userParameterLongType(
            'xe_dissolved_offset_frequency')
        offFreq.value = float(rls.header['sin']['acq_gamma'][0][0]) / H1_GAMMA * float(
            exp.H1resonanceFrequency_Hz) * data_set_config.xe_dissolved_offset_ppm / 1000000
        userParams.userParameterLong.insert(0, offFreq)

    # Remove interleaving
    temp_max_spokes_per_intlv = 0
    if data_set_config.data_type == DataType.DIXON or DataType.UTE:
        temp_max_spokes_per_intlv = enc.encodingLimits.kspace_encoding_step_1.maximum + 1
        enc.encodingLimits.kspace_encoding_step_1.maximum = (
            (enc.encodingLimits.kspace_encoding_step_2.maximum+1) * temp_max_spokes_per_intlv) - 1
        enc.encodingLimits.kspace_encoding_step_2.maximum = 0

    # Save updated header
    header.userParameters = userParams
    header.sequenceParameters = pars
    header.acquisitionSystemInformation = sysInfo
    header.experimentalConditions = exp
    enc.trajectoryDescription = trajDescr
    header.encoding[0] = enc
    header.studyInformation = studyInfo
    header.subjectInformation = subjectInfo

    # account for switching of data labels
    # repetitions are instead contrast with proton/gas/dissolved as 0/1/2
    if data_set_config.data_type != DataType.UTE:
        header.encoding[0].encodingLimits.contrast.minimum = min(
            data_set_config.contrast_order)
        header.encoding[0].encodingLimits.contrast.maximum = max(
            data_set_config.contrast_order)
        header.encoding[0].encodingLimits.contrast.center = 1
        header.encoding[0].encodingLimits.repetition.minimum = 0
        header.encoding[0].encodingLimits.repetition.maximum = 0
        header.encoding[0].encodingLimits.repetition.center = 0
    elif data_set_config.data_type == DataType.UTE:
        header.encoding[0].encodingLimits.contrast.minimum = 0
        header.encoding[0].encodingLimits.contrast.maximum = 0
        header.encoding[0].encodingLimits.contrast.center = 0

    # echoes are no longer contrast but instead sets; contrast limits were updated above
    header.encoding[0].encodingLimits.set.minimum = 1
    header.encoding[0].encodingLimits.set.maximum = int(
        rls.header['sin']['nr_echoes'][0][0])
    header.encoding[0].encodingLimits.set.center = 1

    # calibration doesn't do any encoding
    if data_set_config.data_type == DataType.CALIBRATION:
        header.encoding[0].encodedSpace.fieldOfView_mm.x = np.inf
        header.encoding[0].encodedSpace.fieldOfView_mm.y = np.inf
        header.encoding[0].encodedSpace.fieldOfView_mm.z = np.inf
        header.encoding[0].reconSpace.matrixSize.x = 1
        header.encoding[0].reconSpace.matrixSize.y = 1
        header.encoding[0].reconSpace.matrixSize.z = 1
        header.encoding[0].reconSpace.fieldOfView_mm.x = np.inf
        header.encoding[0].reconSpace.fieldOfView_mm.y = np.inf
        header.encoding[0].reconSpace.fieldOfView_mm.z = np.inf
        header.encoding[0].encodingLimits.kspace_encoding_step_0.maximum = 0
        header.encoding[0].encodingLimits.kspace_encoding_step_0.minimum = 0
        header.encoding[0].encodingLimits.kspace_encoding_step_1.maximum = 0
        header.encoding[0].encodingLimits.kspace_encoding_step_1.minimum = 0

    # add additional encodings for 2nd echo in 1pt dixon
    if data_set_config.multi_echo:
        # all echoes after first will be full projections
        for echo_idx in range(1, int(rls.header['sin']['nr_echoes'][0][0])-1):
            header.encoding.append(copy.deepcopy(header.encoding[0]))
            header.encoding[echo_idx].encodingLimits.kspace_encoding_step_0.minimum = - \
                (header.encoding[echo_idx].encodingLimits.kspace_encoding_step_0.maximum+1)
            header.encoding[echo_idx].encodedSpace.matrixSize.x = int(
                header.encoding[0].encodedSpace.matrixSize.x * 2)

    # finish header update
    dset.write_xml_header(mrd.xsd.ToXML(header))

    # update data headers
    for acqnum in range(dset.number_of_acquisitions()):
        acq_temp = dset.read_acquisition(acqnum)

        # set flag for bonus spectra (prior to reusing set for echoes)
        if data_set_config.bonus_spec == True:
            acq_temp.measurement_uid = acq_temp.idx.set

        if data_set_config.ext_traj == True:
            try:
                if acq_temp.idx.contrast == 0:
                    traj = crds[acq_temp.idx.kspace_encode_step_2,
                                acq_temp.idx.kspace_encode_step_1, :, :]
                else:
                    traj = crds_flyback[acq_temp.idx.kspace_encode_step_2,
                                        acq_temp.idx.kspace_encode_step_1, :, :]
                acq_temp.traj[:] = traj
                acq_temp.sample_time_us = dwell.value
            except:
                pass

        # mrd contrasts = xemrd sets
        acq_temp.idx.set = acq_temp.idx.contrast + 1  # contrast/echo 0 is set 1

        # xemrd contrasts = proton/gas/dissolved 0/1/2
        if data_set_config.data_type == DataType.CALIBRATION:
            # First acqs are dissolved; then gas
            if acqnum < data_set_config.cal_diss_acqs:
                acq_temp.idx.contrast = 2  # 2 = dissolved
            elif acqnum >= data_set_config.cal_diss_acqs:
                acq_temp.idx.contrast = 1  # 1 = gas
            acq_temp.idx.kspace_encode_step_1 = 0
        if data_set_config.data_type == DataType.DIXON:
            acq_temp.idx.contrast = data_set_config.contrast_order[acq_temp.idx.repetition]
        if data_set_config.data_type == DataType.UTE:
            acq_temp.idx.contrast = 0

        acq_temp.idx.repetition = 0  # repetition used for gas vs diss is now unused

        # remove interleaving
        if data_set_config.data_type == DataType.DIXON:
            acq_temp.idx.kspace_encode_step_1 = acq_temp.idx.kspace_encode_step_2 * \
                temp_max_spokes_per_intlv + acq_temp.idx.kspace_encode_step_1
            acq_temp.idx.kspace_encode_step_2 = 0

        # Replace old acq header
        dset.write_acquisition(acq_temp, acqnum)

    # Save dset
    dset.close()

    # Rename
    if data_set_config.data_type == DataType.CALIBRATION:
        try:
            os.remove(os.path.join(mrdName.parent,
                      patientID+'_calibration.h5'))
        except:
            pass
        os.rename(mrdName, os.path.join(
            mrdName.parent, patientID+'_calibration.h5'))
    if data_set_config.data_type == DataType.DIXON:
        try:
            os.remove(os.path.join(mrdName.parent, patientID+'_dixon.h5'))
        except:
            pass
        os.rename(mrdName, os.path.join(mrdName.parent, patientID+'_dixon.h5'))
    if data_set_config.data_type == DataType.UTE:
        try:
            os.remove(os.path.join(mrdName.parent, patientID+'_proton.h5'))
        except:
            pass
        os.rename(mrdName, os.path.join(
            mrdName.parent, patientID+'_proton.h5'))


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="A script/function to convert Xe gas exchange data collected on Philips to the Xe CTC MRD standard.")

    # Define optional arguments
    parser.add_argument(
        "-d", "--data_file",
        type=str,
        default=None,
        help="The path to the data file (default: None)"
    )
    parser.add_argument(
        "-r", "--raw_file",
        type=str,
        default=None,
        help="The path to the raw file (default: None)"
    )
    parser.add_argument(
        "-t", "--traj_file",
        type=str,
        default=None,
        help="The path to the trajectory file (default: None)"
    )

    args = parser.parse_args()

    # Call the main function with parsed arguments
    Gx2XeCTCMRD(data_file=args.data_file,
                raw_file=args.raw_file, traj_file=args.traj_file)
