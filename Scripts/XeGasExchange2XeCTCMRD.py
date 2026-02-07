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
import importlib
from nmr_timefit import NMR_TimeFit
from nmr_mix import NMR_Mix
from scipy.fft import fft, fftshift
import matplotlib.pyplot as plt

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

# Constants
H1_GAMMA = 42577.4688

def make_filtered_mrd_by_indices(src_path, dst_path, keep_indices):
    src = mrd.Dataset(str(src_path), "dataset", create_if_needed=False)
    dst = mrd.Dataset(str(dst_path), "dataset", create_if_needed=True)

    dst.write_xml_header(src.read_xml_header())

    for i in keep_indices:
        dst.append_acquisition(src.read_acquisition(int(i)))

    src.close()
    dst.close()

def find_bonus_spectra_index(acqs, min_ratio=2.0):
    """
    Bonus spectra has a much larger number_of_samples than imaging projections.
    Returns index of the first acquisition that is >= min_ratio * median(others).
    """
    ns = np.array([a.number_of_samples for a in acqs], dtype=float)
    med = np.median(ns)
    # candidates: clearly larger than typical imaging readout
    cand = np.where(ns >= (min_ratio * med))[0]
    return int(cand[0]) if cand.size else None

def movmean(x: np.ndarray, n: int) -> np.ndarray:
    """Compute moving mean of x over n points.

    Args:
        x (np.ndarray): input data of shape (n,)
        n (int): number of points to average over

    Returns:
        np.ndarray: moving mean of shape (n,)
    """
    return np.convolve(x, np.ones((n,)) / n, mode="same")

def gas_phase_contamination_removal(
    data_dissolved: np.ndarray,
    data_gas: np.ndarray,
    sample_time: float,
    freq_gas_acq_diss: float,
    phase_gas_acq_diss: float,
    area_gas_acq_diss: float,
    fa_gas: float,
) -> np.ndarray:
    """Remove gas phase contamination in dissolved k-space.

    Takes gas phase k-space and modifies it using NMR fits and gas phase k0
    to produce the expected gas phase contamination k-space data which is
    then removed from the initial contaminated dissolved phase k-space.

    Args:
        data_dissolved (np.ndarray): dissolved k-space data of shape
            (n_projections, n_points)
        data_gas (np.ndarray): gas phase k-space data of shape
            (n_projections, n_points)
        sample_time (float): dwell time in seconds.
        freq_gas_acq_diss (float): gas frequency offset in dissolved
            spectra acquisition in Hz.
        phase_gas_acq_diss (float): gas phase in dissolved spectra acquisition.
            in degrees.
        area_gas_acq_diss (float): gas area in dissolved spectra acquisition.
        fa_gas (float): gas flip angle in degrees.
    Returns:
        Gas phase corrected dissolved k-space data of shape (n_projections, n_points)
    Author: Matt Willmering
    Paper: https://pubmed.ncbi.nlm.nih.gov/33665905/
    """
    # step 0: calculate parameters
    arr_t = sample_time * np.arange(data_dissolved.shape[1])
    # step 1: modulate contamination (gas) to dissolved frequency - first order
    # phase approximation
    phase_shift1 = 2 * np.pi * freq_gas_acq_diss * arr_t  # calculate phase accumulation
    contamination_kspace1 = data_gas * np.exp(1j * phase_shift1)
    # step 2: zero order phase shift of contamination estimation
    phase_shift2 = phase_gas_acq_diss - 180 / np.pi * np.mean(np.angle(data_gas[:, 0]))
    contamination_kspace2 = contamination_kspace1 * np.exp(
        1j * np.pi / 180 * phase_shift2
    )
    # step 3: scale contamination estimation
    scale_factor = area_gas_acq_diss / movmean(np.abs(data_gas[:, 0]), 100)[-1]
    contamination_kspace3 = (
        contamination_kspace2 * scale_factor / np.cos(np.pi / 180 * fa_gas)
    )
    # step 4: return subtracted contamination
    return data_dissolved - contamination_kspace3

def Gx2XeCTCMRD(data_file=None, raw_file=None, traj_file=None):
    # Get paths
    if data_file == '' and raw_file == '':
        raise RuntimeError(
            f'No raw data file passed. To use file selection dialog, pass None.')
    root = tk.Tk()
    root.withdraw()
    if data_file == None:
        dlName = filedialog.askopenfilename(title='Select .data file', filetypes=[
            ("Philips .data file", "*.data")])
    elif data_file == None:
        dlName = None
    else:
        dlName = data_file
        
    outDir = Path(dlName).parent.absolute()
    if raw_file == None:
        rlsName = filedialog.askopenfilename(title='Select .raw file', filetypes=[
            ("Philips .raw file", "*.raw")], initialdir=outDir)
    elif raw_file == None:
        raw_file = None
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
 
    n_acq_total = dset.number_of_acquisitions()
    n_acq_use = n_acq_total - (1 if data_set_config.exclude_bonus_spec else 0)
    print("n_acq_total =", n_acq_total, "n_acq_use =", n_acq_use,
        "exclude_bonus_spec =", data_set_config.exclude_bonus_spec)

    # If requested, create a new MRD file that excludes the bonus_spec acquisition
    if data_set_config.exclude_bonus_spec:
        mrdPath = Path(mrdName)
        filteredPath = mrdPath.with_name(mrdPath.stem + "_noBonus.h5")

        # build filtered file
        acqs_init      = [None] * n_acq_use
        for acqnum in range(n_acq_use):
            acq_temp = dset.read_acquisition(acqnum)
            acqs_init[acqnum] = acq_temp
        
        bonus_idx = find_bonus_spectra_index(acqs_init)
        print('bonus_idx = ', bonus_idx)
        keep_idx = [i for i in range(dset.number_of_acquisitions()) if i != bonus_idx]
        make_filtered_mrd_by_indices(mrdPath, filteredPath, keep_idx)

        # close original and reopen filtered for the rest of the pipeline
        dset.close()
        mrdName = filteredPath
        dset = mrd.Dataset(str(mrdName), "dataset", create_if_needed=False)

        print("Using filtered MRD:", mrdName, "acqs =", dset.number_of_acquisitions())
    else:
        print("Using original MRD:", mrdName, "acqs =", dset.number_of_acquisitions())


    # Get path to sin file trajectories if necessary
    if data_set_config.ext_traj == True:
        if traj_file == None:
            # Directory containing the running script
            script_dir = Path(__file__).resolve().parent
            repo_root = script_dir.parent

            if 'FLORET'.lower() in rls.header['sin']['scan_name'][0][0].lower():
                ext_coords = repo_root / "resources" / "20251205_083834_Xenon_3D_FLORET_Dixon-NoSpectra.sin"
            else:
                ext_coords = repo_root / "resources" / "20221025_144528_DukeIPF_Gas_Exchange.sin"

            # Safety check
            if not ext_coords.exists():
                ext_coords = filedialog.askopenfilename(title='Select trajectory .sin file', filetypes=[
                    ("Trajectory .sin file", "*.sin")], initialdir=outDir)
        else:
            ext_coords = traj_file
            
        inputDataTraj = rp.PhilipsData(ext_coords)
        inputDataTraj.trajtype = data_set_config.trajorder
        inputDataTraj.delay = data_set_config.gr_delay
        inputDataTraj.readParamOnly = True
        inputDataTraj.compute()
        try:
            traj_type = int(
                inputDataTraj.header['sin']['k_space_traj_type'][0][0])
        except:
            traj_type = 0
        print('traj_type = ', traj_type)

        if traj_type == 1:  # radial
            crds = inputDataTraj.radparams['COORDS']
            if int(inputDataTraj.header['sin']['nr_echoes'][0][0]) > 1:
                crds_flyback = inputDataTraj.radparams['COORDS_FLYBACK']
        elif traj_type == 2:  # spiral
            crds = inputDataTraj.spparams['COORDS_EXPANDED']

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
            inputDataTraj.header['sin']['sample_time_interval'][0][0])
    else:
        dwell.value = float(rls.header['sin']['sample_time_interval'][0][0])
    trajDescr.userParameterDouble.insert(0, dwell)

    ramp_time = mrd.xsd.userParameterLongType('ramp_time')
    ramp_time.value = 0
    if data_set_config.ext_traj == True:
        try:
            ramp_time.value = round(
                float(inputDataTraj.header['sin']['non_cart_fid_slope'][0][0])*dwell.value)
        except:
            pass
    else:
        try:
            ramp_time.value = round(
                float(rls.header['sin']['non_cart_fid_slope'][0][0])*dwell.value)
        except:
            pass
    trajDescr.userParameterLong.insert(0, ramp_time)

    # Determine Gas contamination removal
    if data_set_config.data_type == DataType.DIXON:
        if data_set_config.gas_contam_removal == True:
            if data_set_config.bonus_spec == False:
                print('Contamination removal not possible without bonus spectra - Turning it to false')
                data_set_config.gas_contam_removal = False
        gas_contam_removed = mrd.xsd.userParameterLongType('gas_contam_removed')
        gas_contam_removed.value = int(data_set_config.gas_contam_removal)  # ensure 0/1
        userParams.userParameterLong.insert(0, gas_contam_removed)        

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
        centFreq.value = int(np.round(
            float(rls.header['sin']['acq_gamma'][0][0]) / H1_GAMMA * float(exp.H1resonanceFrequency_Hz)))
        userParams.userParameterLong.insert(0, centFreq)
        offFreq = mrd.xsd.userParameterLongType(
            'xe_dissolved_offset_frequency')
        offFreq.value = int(np.round(float(rls.header['sin']['acq_gamma'][0][0]) / H1_GAMMA * float(
            exp.H1resonanceFrequency_Hz) * data_set_config.xe_dissolved_offset_ppm / 1000000))
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

    # Create data matrices
    Nkx = header.encoding[0].encodingLimits.kspace_encoding_step_0.maximum - header.encoding[0].encodingLimits.kspace_encoding_step_0.minimum + 1
    #acqs = data = contrasts = sets = kzs = kys = kxs = [None] * dset.number_of_acquisitions()
    #acqs = data = contrasts = sets = kzs = kys = kxs = [None] * n_acq_use
    acqs      = [None] * n_acq_use
    data      = [None] * n_acq_use
    contrasts = [None] * n_acq_use
    sets      = [None] * n_acq_use
    kzs       = [None] * n_acq_use
    kys       = [None] * n_acq_use
    kxs       = [None] * n_acq_use

    # Read/store all data/labels for ease/simplicity
    for acqnum in range(n_acq_use):
        acq_temp = dset.read_acquisition(acqnum)
        acqs[acqnum] = acq_temp
        data[acqnum] = acq_temp.data
        contrasts[acqnum] = acq_temp.idx.contrast
        sets[acqnum] = acq_temp.idx.set
        kzs[acqnum] = acq_temp.idx.kspace_encode_step_2
        kys[acqnum] = acq_temp.idx.kspace_encode_step_1
        kxs[acqnum] = acq_temp.number_of_samples

    '''  
    # Get info necessary for gas contamination removal
    if data_set_config.gas_contam_removal == True:
        return
        # see if last projection
        # see if first 
        # see if 2nd mix (scale correction)
        # calculate beta [(last gas phase proj k0 / gas phase area in bonus spec) * (1 / cos(gas flip angle))]
        # calculate dTheta [gas phase in spectra - gas phase in image] 
        # calcute readout times and frequency offset
      # caculate full readout scaling
    ''' 
    # ---------------- Gas contamination removal ----------------
    bonus_idx = find_bonus_spectra_index(acqs)
    if bonus_idx is None:
        print("No bonus spectra found by size; skipping gas contamination removal.")
        data_set_config.gas_contam_removal = False
    else:
        bonus_acq = acqs[bonus_idx]    
    if data_set_config.gas_contam_removal:
        # a) identify bonus spectra acquisition 
        bonus_fid = np.squeeze(bonus_acq.data) 
        bonus_fid = bonus_fid.reshape(-1)     
        print('bonus_fid size =', bonus_fid.size)
        dwell_s = float(dwell.value) * 1e-3 
        t = np.arange(bonus_fid.size) * dwell_s

        # b) fit dissolved spectra 
        plot_fit = True
        disData1_avg = bonus_fid
        fitObj = NMR_TimeFit(
            ydata=disData1_avg,
            tdata=t,
            area=np.array([1, 1, 1]),
            freq=np.array([0, -720, -7700]),
            fwhmL=np.array([250, 200, 30]),
            fwhmG=np.array([0, 200, 0]),
            phase=np.array([0, 0, 0]),
            line_broadening=0,
            zeropad_size=np.size(t),
            method="voigt"
        )
        disfitObj = fitObj.calc_time_fit_residual((-np.inf, np.inf)).T  # [nComp, nParams]

        if plot_fit:
            # --- Measured dissolved spectrum (magnitude) ---
            DisspectralDomainSignal = np.abs(dwell_s * fftshift(fft(disData1_avg)))

            # Frequency axis (Hz)
            N = disData1_avg.size
            f = np.linspace(-0.5, 0.5, N + 1) / dwell_s
            f = f[:-1]

            # --- Build fitted dissolved spectrum from disfitObj ---
            # disfitObj rows: [area, freq, fwhmL, fwhmG, phase] per component
            DisSpectra = NMR_Mix(
                area=disfitObj[:, 0],
                freq=disfitObj[:, 1],
                phase=disfitObj[:, 4],
                fwhmL=disfitObj[:, 2],
                fwhmG=disfitObj[:, 3],
                method="voigt",
            )

            # Time-domain fit (sum components), then FFT -> spectral magnitude
            time = dwell_s * np.arange(N)
            DissFit_td = DisSpectra.calcComponentTimeDomainSignal(time)         # [N, nComp] complex
            DissFit_spec = dwell_s * fftshift(fft(DissFit_td, axis=0), axes=0)  # [N, nComp] complex
            DissFit = np.abs(np.sum(DissFit_spec, axis=1))                      # [N] magnitude

            # --- Plot (same style you showed) ---
            fig, ax3 = plt.subplots()

            # ---- select central 25% of the spectrum ----
            N = len(f)
            center = N // 2
            half_width = int(0.25 * N / 2)   # 25% total → 12.5% on each side

            idx = slice(center - half_width, center + half_width)

            # ---- plot only central portion ----
            ax3.plot(f[idx], np.flip(DisspectralDomainSignal)[idx], 'bo', markerfacecolor='b')
            ax3.plot(f[idx], np.flip(DissFit)[idx], '-r')

            ax3.set_title('Dissolved Phase Spectrum (Central 25%)')
            ax3.set_xlabel('Frequency (Hz)')
            ax3.set_ylabel('NMR Signal Intensity (a.u)')

            plt.tight_layout()
            plt.show()

        gas_idx = 2  # the -7700 Hz component
        GasArea_DissAcq = float(disfitObj[gas_idx, 0])
        GasPhase_DissAcq_deg = float(disfitObj[gas_idx, 4])
        FreqOffset_Hz = float(disfitObj[gas_idx, 1])

        # c) build gas/diss kspace matrices from imaging acquisitions
        # after relabeling later: contrast 1 = gas, 2 = dissolved
        # In Philips MRD: repetition distinguishes gas vs dissolved for Dixon
        gas_rep  = data_set_config.contrast_order.index(1)   # rep that maps to Xe gas (contrast=1)
        diss_rep = data_set_config.contrast_order.index(2)   # rep that maps to Xe dissolved (contrast=2)

        gas_acq_idx  = [i for i,a in enumerate(acqs) if i != bonus_idx and a.idx.repetition == gas_rep]
        diss_acq_idx = [i for i,a in enumerate(acqs) if i != bonus_idx and a.idx.repetition == diss_rep]

        GasKSpace  = np.stack([np.squeeze(acqs[i].data).reshape(-1) for i in gas_acq_idx],  axis=1)   # [RO, Nproj]
        DissKSpace = np.stack([np.squeeze(acqs[i].data).reshape(-1) for i in diss_acq_idx], axis=1)  # [RO, Nproj]

        # d) correct dissolved kspace
        DissKSpace_corr = gas_phase_contamination_removal(
            DissKSpace=DissKSpace,
            GasKSpace=GasKSpace,
            dwell_s=dwell_s,
            FreqOffset_Hz=FreqOffset_Hz,
            GasPhase_DissAcq_deg=GasPhase_DissAcq_deg,
            GasArea_DissAcq=GasArea_DissAcq,
            GasFA_deg=float(data_set_config.flip_angle_gas),
        )

        # e) write corrected data back into dissolved acquisitions
        for col, i in enumerate(diss_acq_idx):
            acqs[i].data[:] = DissKSpace_corr[:, col].reshape(acqs[i].data.shape)
            
        print("Gas contamination removal applied.")
    # -----------------------------------------------------------

    # update data and headers
    # Version   |   
    # Philips:  | location | average | extr1   | mix      | card  | dynamic    | echo     | kz | ky | kx
    # MRD:      | slice    | average | segment | set      | phase | repetition | contrast | kz | ky | kx
    # XeCTCMRD: | N/A      | N/A     | N/A     | spec/img | N/A   | contrast   | set      | kz | ky | kx
    for acqnum in range(len(acqs)):
        acq_temp = acqs[acqnum]

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
        
    if data_set_config.gas_contam_removal:
        # remove bonus spectra acquisition from the final file (exclude bonus_idx)
        keep_idx = [i for i in range(dset.number_of_acquisitions()) if i != bonus_idx]
        mrdPath = Path(mrdName)
        filteredPath = Path(rlsName).with_suffix('').with_name(Path(rlsName).stem + "_noBonus.h5")

        dset.close()  # close before copying
        make_filtered_mrd_by_indices(mrdPath, filteredPath, keep_idx)

        mrdName = filteredPath
        dset = mrd.Dataset(str(mrdName), "dataset", create_if_needed=False)
        print("Using filtered MRD:", mrdName, "acqs =", dset.number_of_acquisitions())

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
            if data_set_config.gas_contam_removal:
                os.remove(os.path.join(mrdName.parent, patientID+'_dixon_corrected.h5'))
            else:
                os.remove(os.path.join(mrdName.parent, patientID+'_dixon.h5'))
        except:
            pass
        if data_set_config.gas_contam_removal:
            os.rename(mrdName, os.path.join(mrdName.parent, patientID+'_dixon_corrected.h5'))
        else:
            os.rename(mrdName, os.path.join(mrdName.parent, patientID+'_dixon.h5'))
    if data_set_config.data_type == DataType.UTE:
        try:
            os.remove(os.path.join(mrdName.parent, patientID+'_proton.h5'))
        except:
            pass
        os.rename(mrdName, os.path.join(
            mrdName.parent, patientID+'_proton.h5'))

''' 
if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="A script/function to convert Xe gas exchange data collected on Philips to the Xe CTC MRD standard.")

    # Define optional arguments
    parser.add_argument(
        "-d", "--data_file",
        type=str,
        default=None,
        help="The path to the data file (default: None). Pass '' if want file dialog."
    )
    parser.add_argument(
        "-r", "--raw_file",
        type=str,
        default=None,
        help="The path to the raw file (default: None). Pass '' if want file dialog."
    )
    parser.add_argument(
        "-t", "--traj_file",
        type=str,
        default=None,
        help="The path to the trajectory file (default: None). Pass '' if want file dialog."
    )

    args, unknown_args = parser.parse_known_args()

    # Call the main function with parsed arguments
    Gx2XeCTCMRD(data_file=args.data_file,
                raw_file=args.raw_file, traj_file=args.traj_file)
'''
    
    
data_file = r"D:\OneDrive - cchmc\Lab\Random Subject analysis\Philips2MRD_data\UF_3D_radial_Dixon\6)Gas_Exchange\raw_409.data"
raw_file = r"D:\OneDrive - cchmc\Lab\Random Subject analysis\Philips2MRD_data\UF_3D_radial_Dixon\6)Gas_Exchange\20251124_162012_Xenon_3D_Radial_Dixon.raw"
traj_file = None 
Gx2XeCTCMRD(data_file, raw_file, traj_file)
