
#%Import libraries
import sys  # noqa: E402
from pathlib import Path  # noqa: E402
p2mDir = Path(__file__).parent.parent.absolute()  # noqa: E402
sys.path.append(str(p2mDir))  # noqa: E402
import philips2mrd as p2m
# Import read philips script (non-compiled)
import readphilips.ReadPhilips as rp
import ismrmrd as mrd
from Scripts.XeGasExchange2XeCTCMRD_Config import Config, DataType
import argparse
from tkinter import filedialog
import tkinter as tk
import numpy as np
import copy
import os
import math
import importlib
from nmr_timefit import NMR_TimeFit
from nmr_mix import NMR_Mix
from scipy.fft import fft, fftshift
import matplotlib.pyplot as plt


if sys.version_info.major != 3:
    raise RuntimeError('Requires python 3')

#Define Functions

major_version = sys.version_info.major
minor_version = sys.version_info.minor
rp_name = f"rp.rp{major_version}{minor_version}"
# try:
#     rp = importlib.import_module(rp_name)
# except ModuleNotFoundError:
#     raise RuntimeError(
#         f'ReadPhilips not compiled for Python {sys.version_info.major}.{sys.version_info.minor}')

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

def find_odd_spectra_indices(acqs, min_ratio=2.0):
    """
    Identify 'odd' acquisitions by number_of_samples (vector size).

    Strategy
    - Compute number_of_samples for each acquisition.
    - Split into the dominant (larger count) group and the minority group based on frequency.
      (Assumption: imaging readouts dominate; bonus spectra are few and have larger size.)
    - Use the dominant group's median as the reference.
    - Flag odds as those that deviate strongly from that median:
        - large odds: >= min_ratio * median(dominant)
        - small odds: <= median(dominant) / min_ratio  (kept for completeness)

    Returns
    - odd_indices: list[int] of all odd acquisitions
    - info: dict with counts/medians for debugging
    """
    ns = np.array([int(a.number_of_samples) for a in acqs], dtype=int)

    # 1) Find unique sizes + counts
    sizes, counts = np.unique(ns, return_counts=True)
    if sizes.size <= 1:
        return [], {
            "unique_sizes": sizes.tolist(),
            "counts": counts.tolist(),
            "median_ref": float(np.median(ns)),
            "n_odd": 0
        }

    # 2) Dominant size group = the most frequent size (assumed "normal")
    dominant_size = int(sizes[np.argmax(counts)])
    dominant_mask = (ns == dominant_size)

    # 3) Reference median from the dominant group (robust, and uses the majority)
    ref_med = float(np.median(ns[dominant_mask]))

    # 4) Odd candidates: anything that differs in size AND is far from ref median
    different_size = (ns != dominant_size)
    odd_large = different_size & (ns >= int(np.ceil(min_ratio * ref_med)))
    odd_small = different_size & (ns <= int(np.floor(ref_med / min_ratio)))

    odd_mask = odd_large | odd_small
    odd_indices = np.where(odd_mask)[0].tolist()

    info = {
        "unique_sizes": sizes.tolist(),
        "counts": counts.tolist(),
        "dominant_size": dominant_size,
        "n_dominant": int(dominant_mask.sum()),
        "n_different_size": int(different_size.sum()),
        "median_ref": ref_med,
        "n_odd": int(len(odd_indices)),
        "odd_indices": odd_indices,
        "odd_sizes": ns[odd_mask].tolist(),
    }
    return odd_indices, info

def plot_complex_spectra_panels(data, fs=None, do_fft=False, suptitle='Spectra'):
    """
    data : complex array (n_spectra, N) or (N, n_spectra)
    fs   : sampling frequency (Hz) if FFT desired
    do_fft : apply FFT before plotting
    """

    x = np.asarray(data)

    # Ensure shape = (n_spectra, N)
    if x.ndim != 2:
        raise ValueError("Input must be 2D array")

    if x.shape[0] < x.shape[1]:
        n_spec = x.shape[0]
    else:
        x = x.T
        n_spec = x.shape[0]

    # FFT if needed
    if do_fft:
        x = np.fft.fftshift(np.fft.fft(x, axis=1), axes=1)

    N = x.shape[1]

    # X axis
    if do_fft and fs is not None:
        xx = np.fft.fftshift(np.fft.fftfreq(N, d=1/fs))
        xlabel = 'Frequency (Hz)'
    else:
        xx = np.arange(N)
        xlabel = 'Sample'

    # Determine subplot grid
    ncols = min(3, n_spec)
    nrows = math.ceil(n_spec / ncols)

    fig, axes = plt.subplots(nrows, ncols, figsize=(5*ncols, 3*nrows), squeeze=False)

    for i in range(n_spec):
        r = i // ncols
        c = i % ncols

        ax = axes[r, c]

        ax.plot(xx, np.abs(x[i]), label='Magnitude')
        ax.set_title(f"{'Pre' if i < n_spec/2 else 'Post'} {i+1}")
        ax.set_xlabel(xlabel)
        ax.set_ylabel('|S|')
        ax.grid(alpha=0.3)

    # Remove unused panels
    for j in range(i+1, nrows*ncols):
        r = j // ncols
        c = j % ncols
        fig.delaxes(axes[r][c])

    fig.suptitle(suptitle)
    plt.tight_layout()
    plt.show()

def movmean2(x: np.ndarray, n: int) -> np.ndarray:
    ''' 
    Compute moving mean of x over n points.
    Args:
        x (np.ndarray): input data of shape (n,)
        n (int): number of points to average over
    Returns:
        np.ndarray: moving mean of shape (n,)
    '''
    return np.convolve(x, np.ones((n,)) / n, mode="same")

def movmean(x: np.ndarray, n: int, axis: int = 0) -> np.ndarray:
    """
    Moving mean along a given axis, MATLAB-like 'movmean(..., n)'.
    Uses convolution with 'same' length along that axis.
    Works for complex arrays and N-D arrays.
    """
    x = np.asarray(x)
    if n <= 1:
        return x

    kernel = np.ones(n, dtype=np.float64) / n

    # Apply 1D convolution along the specified axis
    return np.apply_along_axis(lambda v: np.convolve(v, kernel, mode="same"), axis, x)

def downsample(x: np.ndarray, factor: int, axis: int = 0) -> np.ndarray:
    """
    Downsample by integer factor along a specified axis (MATLAB downsample behavior).
    """
    if factor <= 1:
        return x
    slc = [slice(None)] * x.ndim
    slc[axis] = slice(None, None, factor)
    return x[tuple(slc)]

def gas_phase_contamination_removal(
    data_dissolved: np.ndarray,
    data_gas: np.ndarray,
    sample_time: float,
    freq_gas_acq_diss: float,
    phase_gas_acq_diss: float,
    area_gas_acq_diss: float,
    fa_gas: float,
) -> np.ndarray:
    '''
    Remove gas phase contamination in dissolved k-space.
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
    ''' 
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

def make_timefit_x0_feasible(obj, lb, ub):
    """
    Modifies obj.area/freq/fwhmL/fwhmG/phase so that the flattened x0
    lies inside [lb, ub]. Uses small eps offsets to avoid equality issues.
    """
    eps = 1e-12

    x0 = np.array([obj.area, obj.freq, obj.fwhmL, obj.fwhmG, obj.phase]).flatten()

    # If you used epsilon bounds for "fixed" params, keep x0 below ub-eps
    x0 = np.minimum(np.maximum(x0, lb), ub)
    # Nudge off exact upper bound just in case (helps with inf not needed)
    finite_ub = np.isfinite(ub)
    x0[finite_ub] = np.minimum(x0[finite_ub], ub[finite_ub] - eps)

    # Write back into the object (order must match x0 construction)
    x0m = x0.reshape(5, -1)   # [area; freq; fwhmL; fwhmG; phase] each length 3
    obj.area  = x0m[0, :]
    obj.freq  = x0m[1, :]
    obj.fwhmL = x0m[2, :]
    obj.fwhmG = x0m[3, :]
    obj.phase = x0m[4, :]

    return obj

def plot_traj(crds, n_traj=10):
    """
    crds shape assumed (..., 3)  e.g. (10,95,58,3)
    n_traj: number of trajectories to display
    """
    traj = crds.reshape(-1, crds.shape[-2], 3)   # (Ntraj, Nsamp, 3)
    n_traj = min(n_traj, len(traj))

    ax = plt.figure().add_subplot(111, projection="3d")
    colors = plt.cm.jet(np.linspace(0,1,n_traj))

    for i in range(n_traj):
        ax.plot(*traj[i].T, color=colors[i], lw=1)

    plt.show()

#%Import data
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

    # Run Config
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

    acqs_init      = [None] * n_acq_total
    for acqnum in range(n_acq_total):
        acq_temp = dset.read_acquisition(acqnum)
        acqs_init[acqnum] = acq_temp
    bonus_idx, bonus_info = find_odd_spectra_indices(acqs_init)
    print('bonus_idx = ', bonus_idx)

    n_acq_use = n_acq_total - (len(bonus_idx) if data_set_config.exclude_bonus_spec else 0)
    print("n_acq_total =", n_acq_total, "n_acq_use =", n_acq_use,
        "exclude_bonus_spec =", data_set_config.exclude_bonus_spec)

    #Create a new MRD file that excludes the bonus_spec acquisition
    if data_set_config.exclude_bonus_spec:

        mrdPath = Path(mrdName)
        filteredPath = mrdPath.with_name(mrdPath.stem + "_noBonus.h5")

        # build filtered file
        keep_idx = [i for i in range(dset.number_of_acquisitions()) 
                    if i not in bonus_idx]
        make_filtered_mrd_by_indices(mrdPath, filteredPath, keep_idx)

        # close original and reopen filtered for the rest of the pipeline
        dset.close()
        mrdName = filteredPath
        dset = mrd.Dataset(str(mrdName), "dataset", create_if_needed=False)
        data_set_config.gas_contam_removal = False # force it here just in case
        print("Using filtered MRD:", mrdName, "acqs =", dset.number_of_acquisitions())
    else:
        print("Using original MRD:", mrdName, "acqs =", dset.number_of_acquisitions())

    # Get path to fixed sin file trajectories if necessary
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

        if traj_type == 1:  # radial
            crds = inputDataTraj.radparams['COORDS']
            if int(inputDataTraj.header['sin']['nr_echoes'][0][0]) > 1:
                crds_flyback = inputDataTraj.radparams['COORDS_FLYBACK']
        elif traj_type == 2:  # spiral
            crds = inputDataTraj.spparams['COORDS_EXPANDED']
    plotting = data_set_config.plotting
    if plotting:
        plot_traj(crds, n_traj=100)

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
        # which method   
        gas_contam_method = mrd.xsd.userParameterStringType('gas_contam_method')
        gas_contam_method.value = data_set_config.gas_contam_method  # ensure 0/1
        userParams.userParameterString.insert(0, gas_contam_method) 

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

    # ---------------- Gas contamination removal ----------------
    data_set_config.gas_contam_removal = False
    if data_set_config.gas_contam_removal:
        bonus_idx, bonus_info = find_odd_spectra_indices(acqs)
        print('bonus_idx = ', bonus_idx)
        print('bonus_info = ', bonus_info)
        if bonus_idx is None:
            print("No bonus spectra found by size; skipping gas contamination removal.")
            data_set_config.gas_contam_removal = False
        else:
            bonus_acqs = [acqs[int(i)] for i in bonus_idx]         
        # a) identify bonus spectra acquisition 
        bonus_fids = [np.squeeze(acq.data) for acq in bonus_acqs]
        # Convert list to numpy array
        bonus_fid = np.array(bonus_fids)
        plotting = True
        if plotting:
            plot_complex_spectra_panels(bonus_fid, do_fft=False, suptitle='Bonus Spectra')

        print('bonus_fid size =', bonus_fid[0].size)
        dwell_s = float(dwell.value) * 1e-6 
        
        PreDissolvedFID = bonus_fid[0]
        PostDissolvedFID = bonus_fid[3]
        PreGasFID = bonus_fid[1]
        PostGasFID = bonus_fid[4]

        # ---- K-space ----
        gas_rep  = data_set_config.contrast_order.index(1)   # repetition that maps to Xe gas (contrast=1)
        diss_rep = data_set_config.contrast_order.index(2)   # repetition that maps to Xe dissolved (contrast=2)

        # bonus_idx is a list -> use a set for O(1) membership
        bonus_set = set(int(i) for i in bonus_idx)

        gas_acq_idx  = [i for i, a in enumerate(acqs)
                        if (i not in bonus_set) and (a.idx.repetition == gas_rep)]

        diss_acq_idx = [i for i, a in enumerate(acqs)
                        if (i not in bonus_set) and (a.idx.repetition == diss_rep)]

        GasKSpaceInit  = np.stack([np.squeeze(acqs[i].data).reshape(-1) for i in gas_acq_idx], axis=1)   # [RO, Nproj]
        DissolvedKSpaceInit = np.stack([np.squeeze(acqs[i].data).reshape(-1) for i in diss_acq_idx], axis=1)  # [RO, Nproj]

        OvsFactor = 4
        extraOvs = True
        if extraOvs:
            # ---- Dissolved FIDs ----
            PreDissolvedFID  = movmean(PreDissolvedFID, OvsFactor)
            PreDissolvedFID  = PreDissolvedFID[::OvsFactor]

            PostDissolvedFID = movmean(PostDissolvedFID, OvsFactor)
            PostDissolvedFID = PostDissolvedFID[::OvsFactor]

            # ---- Gas FIDs ----
            PreGasFID  = movmean(PreGasFID, OvsFactor)
            PreGasFID  = PreGasFID[::OvsFactor]

            PostGasFID = movmean(PostGasFID, OvsFactor)
            PostGasFID = PostGasFID[::OvsFactor]
            
            # ---- K-space ----
            DissolvedKSpaceInit = movmean(DissolvedKSpaceInit, OvsFactor)
            DissolvedKSpaceInit = DissolvedKSpaceInit[::OvsFactor]

            GasKSpaceInit = movmean(GasKSpaceInit, OvsFactor)
            GasKSpaceInit = GasKSpaceInit[::OvsFactor]

            # ---- Adjust dwell time ----
            #dwell_s = dwell_s * OvsFactor            

        if plotting:
            # Ensure at least 3 projections exist
            n_proj_gas = min(3, GasKSpaceInit.shape[1])
            n_proj_diss = min(3, DissolvedKSpaceInit.shape[1])

            x_gas = np.arange(GasKSpaceInit.shape[0])
            x_diss = np.arange(DissolvedKSpaceInit.shape[0])

            fig, axes = plt.subplots(2, 3, figsize=(15, 6), sharex=False)

            for i in range(n_proj_gas):
                axes[0, i].plot(x_gas, np.abs(GasKSpaceInit[:, i]), linewidth=1.5)
                axes[0, i].set_title(f'Gas Projection {i+1}')
                axes[0, i].set_ylabel('|Signal|')
                axes[0, i].grid(alpha=0.3)

            for i in range(n_proj_diss):
                axes[1, i].plot(x_diss, np.abs(DissolvedKSpaceInit[:, i]), linewidth=1.5)
                axes[1, i].set_title(f'Dissolved Projection {i+1}')
                axes[1, i].set_ylabel('|Signal|')
                axes[1, i].set_xlabel('Readout Index')
                axes[1, i].grid(alpha=0.3)

            plt.tight_layout()
            plt.show()
            
            DissolvedKSpaceInit = movmean(DissolvedKSpaceInit, OvsFactor)
            DissolvedKSpaceInit = DissolvedKSpaceInit[::OvsFactor]

            GasKSpaceInit = movmean(GasKSpaceInit, OvsFactor)
            GasKSpaceInit = GasKSpaceInit[::OvsFactor]

            # ---- Adjust dwell time ----
            #dwell_s = dwell_s * OvsFactor

            # ---- Attenuation Correction ----
            # MATLAB: PostGasFID(1,1) and GasKSpaceInit(1,end,end)
        
        scaleFac = np.abs(PostGasFID[0]) / np.abs(GasKSpaceInit[0, -1])
        scaleFac1dis = np.abs(PostDissolvedFID[0]) / np.abs(DissolvedKSpaceInit[0, -1])

        DissolvedKSpaceInit = DissolvedKSpaceInit * scaleFac
        GasKSpaceInit = GasKSpaceInit * scaleFac

        time = dwell_s * np.arange(len(PostDissolvedFID), dtype=float)

        # -------------------- Initial guesses --------------------
        A0 = np.abs(PostDissolvedFID[0])  # MATLAB PostDissolvedFID(1)

        area_guess  = np.array([0.27*A0, 1.0*A0, 0.10*A0], dtype=float)
        freq_guess  = np.array([534.0, -126.0, -7084.0], dtype=float)

        fwhmL_guess = np.array([300.0, 273.0, 44.0], dtype=float)
        fwhmG_guess = np.array([0.0, 275.0, 0.0], dtype=float)

        phase_guess = np.array([-95.0, 151.0, 9.0], dtype=float)  # degrees (your code uses degrees)

        # -------------------- Bounds  --------------------
        area_lb  = np.array([0.0, 0.0, 0.0], dtype=float)
        area_ub  = np.array([1e10, 1e10, 1e10], dtype=float)

        freq_lb  = np.array([500.0, -2000.0, -12000.0], dtype=float)
        freq_ub  = np.array([2000.0, 0.0, -4000.0], dtype=float)

        fwhmL_lb = np.array([0.0, 0.0, 0.0], dtype=float)
        fwhmL_ub = np.array([np.inf, np.inf, np.inf], dtype=float)

        eps = 1e-12  # small number so SciPy sees lb < ub, but parameter is essentially fixed

        fwhmG_lb = np.array([0.0, 0.0, 0.0], dtype=float)
        fwhmG_ub = np.array([eps, np.inf, eps], dtype=float)  # instead of [0, inf, 0]

        phase_lb = np.array([-np.inf, -np.inf, -np.inf], dtype=float)
        phase_ub = np.array([ np.inf,  np.inf,  np.inf], dtype=float)

        lb = np.concatenate([area_lb, freq_lb, fwhmL_lb, fwhmG_lb, phase_lb])
        ub = np.concatenate([area_ub, freq_ub, fwhmL_ub, fwhmG_ub, phase_ub])

        # ensure feasibility before fitting
        PrependedDissolvedNMRFit = make_timefit_x0_feasible(PrependedDissolvedNMRFit, lb, ub)
        PrependedDissolvedNMRFit.fit_time_signal_residual(bounds=(lb, ub))
        
        # -------------------- 2) Fit appended dissolved (Post), seeded by Pre-fit --------------------
        AppendedDissolvedNMRFit = NMR_TimeFit(
            ydata=PostDissolvedFID,
            tdata=time,
            area=area_guess,
            freq=PrependedDissolvedNMRFit.freq,
            fwhmL=PrependedDissolvedNMRFit.fwhmL,
            fwhmG=PrependedDissolvedNMRFit.fwhmG,
            phase=PrependedDissolvedNMRFit.phase,
            line_broadening=0,
            zeropad_size=time.size,
            method="voigt",
        )

        # Optional: initial fit (unbounded) like MATLAB
        AppendedDissolvedNMRFit.fit_time_signal_residual(bounds=(-np.inf, np.inf))

        # Critical: make x0 feasible for bounded refit (SciPy requirement)
        AppendedDissolvedNMRFit = make_timefit_x0_feasible(AppendedDissolvedNMRFit, lb, ub)

        # Bounded refit (MATLAB "setBounds + refit")
        AppendedDissolvedNMRFit.fit_time_signal_residual(bounds=(lb, ub))

        # -------------------- 3) MATLAB: dwell_s * fftshift(fft(calcComponentTimeDomainSignal(time),[],1),1) --------------------
        # Use ONLY calcComponentTimeDomainSignal, then sum components:
        comp_td = AppendedDissolvedNMRFit.calcComponentTimeDomainSignal(time)  # (Nt, 3)
        fit_td  = np.sum(comp_td, axis=1)                                     # (Nt,)

        AppendedDissolvedFit = dwell_s * np.fft.fftshift(np.fft.fft(fit_td))
        phase = AppendedDissolvedNMRFit.phase[2]
        area  = AppendedDissolvedNMRFit.area[2]

        if plotting:
            f = AppendedDissolvedNMRFit.f
            spec = AppendedDissolvedNMRFit.spectral_signal  # measured spectrum

            # If AppendedDissolvedFit is component-resolved (Nfreq, 3)
            fit_components = AppendedDissolvedFit
            fit_sum = np.sum(fit_components, axis=1)

            # ---- Figure ----
            plt.figure(figsize=(16, 9), facecolor="white")

            # ---------------- Data ----------------
            plt.plot(f, np.abs(spec), 'b', label='Spectrum - Magnitude')
            plt.plot(f, np.real(spec), 'r', label='Spectrum - Real')
            plt.plot(f, np.imag(spec), color=(0, 0.6, 0.2), label='Spectrum - Imaginary')

            # ---------------- Fits ----------------
            plt.plot(f, np.abs(fit_sum), color=(0, 0, 1, 0.33), linewidth=3, label='Fit - Magnitude')
            plt.plot(f, np.real(fit_sum), color=(1, 0, 0, 0.33), linewidth=3, label='Fit - Real')
            plt.plot(f, np.imag(fit_sum), color=(0, 0.6, 0.2, 0.33), linewidth=3, label='Fit - Imaginary')

            # ---------------- Components (Real only) ----------------
            plt.fill_between(f, 0, np.real(fit_components[:, 0]),
                            color='r', alpha=0.33, label='RBC - Real')

            plt.fill_between(f, 0, np.real(fit_components[:, 1]),
                            color='b', alpha=0.33, label='Barrier - Real')

            plt.fill_between(f, 0, np.real(fit_components[:, 2]),
                            color=(0, 0.6, 0.2), alpha=0.33, label='Gas - Real')

            # ---------------- Settings ----------------
            plt.legend(loc='best', ncol=3, frameon=False)
            plt.xlim([-8000, 2000])
            plt.gca().invert_xaxis()  # MATLAB 'XDir','reverse'
            plt.xticks(fontsize=18)
            plt.yticks(fontsize=18)

            ratio = AppendedDissolvedNMRFit.area[0] / AppendedDissolvedNMRFit.area[1]

            plt.title(f'Dissolved Phase Spectra and Fit: RBC/Barrier = {ratio:.3f}', fontsize=22)
            plt.xlabel('Frequency (Hz)', fontsize=20)
            plt.ylabel('NMR Signal (a.u.)', fontsize=20)

            plt.tight_layout()
            plt.show()

            print('Fitting Spectrum Completed.')

        freq_jump = 7143
        # d) correct dissolved kspace
        DissKSpace_corr = gas_phase_contamination_removal(
            data_dissolved=DissolvedKSpaceInit,
            data_gas=GasKSpaceInit,
            sample_time=dwell_s,
            freq_gas_acq_diss = -freq_jump,
            phase_gas_acq_diss = phase,
            area_gas_acq_diss = area,
            fa_gas=float(data_set_config.flip_angle_gas),
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
                #acq_temp.traj[:] = traj
                n = min(acq_temp.traj.shape[0], traj.shape[0])
                acq_temp.traj[:n, :] = traj[:n, :].astype(acq_temp.traj.dtype, copy=False)
                                
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

#%
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
    
##Gas data Files
# data_file = r"C:\Users\HUSDQ4\OneDrive - cchmc\cincy_work\all_projects_data_work\gex_analysis\healthy_all\CPIR_protocol\20230316-ILD-HC-066\raw_405.data"
# raw_file = r"C:\Users\HUSDQ4\OneDrive - cchmc\cincy_work\all_projects_data_work\gex_analysis\healthy_all\CPIR_protocol\20230316-ILD-HC-066\20230316_125922_CPIR_Gas_Exchange.raw"
# traj_file = r"C:\Users\HUSDQ4\Desktop\IRC186-507_GX_test\gas\rls_fixed\20200210_133229_Dissolved_Xe_20191008 - 3T-T1.sin"

# fixed_traj_file = r"C:\Users\HUSDQ4\Desktop\IRC186-507_GX_test\gas\rls_fixed\20260402_CPIR_Gas_Exchange_3T-T2.sin"

##Proton/mask data Files
data_file = r"C:\Users\bda5ik\Downloads\raw_405.data"
raw_file = r"C:\Users\bda5ik\Downloads\20240917_120146_CPIR_Gas_Exchange.raw"
traj_file =  r"C:\Users\bda5ik\Downloads\20200210_133229_Dissolved_Xe_20191008 - 3T-T1.sin"



Gx2XeCTCMRD(data_file, raw_file, traj_file)

#%%
