
# main script
import sys
import os
from pathlib import Path
repo_root = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(repo_root))
import philips2mrd as p2m
import readphilips.ReadPhilips as rp
import ismrmrd as mrd
from Scripts.XeGasExchange2XeCTCMRD_Config import Config, DataType
import argparse
from tkinter import filedialog
import tkinter as tk
import numpy as np
import copy
import math
from Scripts.nmr_timefit import NMR_TimeFit
from scipy.fft import fft, fftshift
import matplotlib.pyplot as plt
from scipy.io import savemat

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

def movmean_truncate(x: np.ndarray, n: int, axis: int = 0) -> np.ndarray:
    """
    MATLAB-like moving mean along a given axis.
    Near the boundaries, the window is truncated rather than zero-padded.
    """
    x = np.asarray(x)
    if n <= 1:
        return x

    def movmean_1d(v: np.ndarray, n: int) -> np.ndarray:
        v = np.asarray(v)
        out = np.empty_like(v, dtype=np.result_type(v, np.float64))

        half_left = (n - 1) // 2
        half_right = n // 2

        for i in range(v.size):
            lo = max(0, i - half_left)
            hi = min(v.size, i + half_right + 1)
            out[i] = np.mean(v[lo:hi])

        return out

    return np.apply_along_axis(lambda v: movmean_1d(v, n), axis, x)

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
    # Make sure inputs are arrays
    data_dissolved = np.asarray(data_dissolved)
    data_gas = np.asarray(data_gas)

    # MATLAB expects [RO, Projections]
    if data_dissolved.ndim != 2 or data_gas.ndim != 2:
        raise ValueError("data_dissolved and data_gas must be 2D arrays with shape [RO, Projections].")

    if data_dissolved.shape != data_gas.shape:
        raise ValueError(
            f"data_dissolved and data_gas must have the same shape. "
            f"Got {data_dissolved.shape} and {data_gas.shape}."
        )

    # ------------------------------------------------------------------
    # Step 0: calculate parameters
    # MATLAB: time = dwell_s*(0:size(GasKSpace,1)-1);
    # ------------------------------------------------------------------
    time = sample_time * np.arange(data_gas.shape[0], dtype=float)  # [RO]

    # ------------------------------------------------------------------
    # Step 1: modulate contamination (gas) to dissolved frequency
    # MATLAB:
    # Step1PhaseShift = time' * FreqOffset * 2 * pi;
    # Step1ContamKSpace = GasKSpace .* exp(1i*Step1PhaseShift);
    # ------------------------------------------------------------------
    step1_phase_shift = time[:, None] *freq_gas_acq_diss * 2.0 * np.pi  # [RO, 1]
    step1_contam_kspace = data_gas * np.exp(1j * step1_phase_shift)       # [RO, Projections]

    # ------------------------------------------------------------------
    # Step 2: zero-order phase shift
    # MATLAB:
    # GasPhaseMean = movmean(rad2deg(angle(GasKSpace(1,:))),100);
    # Step2PhaseShift = GasPhase_DissAcq - GasPhaseMean(1,end);
    # Step2ContamKSpace = Step1ContamKSpace * exp(1i*deg2rad(Step2PhaseShift));
    # ------------------------------------------------------------------
    gas_phase_mean = movmean_truncate(np.rad2deg(np.angle(data_gas[0, :])), 100)
    step2_phase_shift = phase_gas_acq_diss - gas_phase_mean[-1]
    step2_contam_kspace = step1_contam_kspace * np.exp(1j * np.deg2rad(step2_phase_shift))

    # ------------------------------------------------------------------
    # Step 3: scale contamination
    # MATLAB:
    # GasK0Mean = movmean(abs(GasKSpace(1,:)),100);
    # GasArea = GasK0Mean(1,end);
    # Step3Scale = GasArea_DissAcq/GasArea;
    # Step3ContamKSpace = Step2ContamKSpace * Step3Scale / cosd(GasFA);
    # ------------------------------------------------------------------
    gas_k0_mean = movmean_truncate(np.abs(data_gas[0, :]), 100)
    gas_area = gas_k0_mean[-1]
    step3_scale = area_gas_acq_diss / gas_area
    step3_contam_kspace = step2_contam_kspace * step3_scale / np.cos(np.deg2rad(fa_gas))

    # ------------------------------------------------------------------
    # Subtract contamination
    # MATLAB:
    # CorrectedDissKSpace = DissKSpace - Step3ContamKSpace;
    # ------------------------------------------------------------------
    corrected_diss_kspace = data_dissolved - step3_contam_kspace

    return corrected_diss_kspace

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

def write_kspace_to_acqs(acqs, acq_idx, kspace_mat, mode="pad"):
    """
    acq_idx: list of acquisition indices
    kspace_mat: array of shape [Nsamp, Nproj]
    mode: 'pad' or 'trim'
    """
    trim_lengths = {}

    # step 1: write/pad into existing acquisition buffers
    for col, i in enumerate(acq_idx):
        src = np.asarray(kspace_mat[:, col]).reshape(-1)
        n = src.size

        nchan = acqs[i].data.shape[0]
        old_nsamp = acqs[i].data.shape[1]

        if n > old_nsamp:
            raise ValueError(
                f"Corrected data longer than acquisition buffer: n={n}, old_nsamp={old_nsamp}, acq={i}"
            )

        tmp = np.zeros((nchan, old_nsamp), dtype=acqs[i].data.dtype)
        tmp[0, :n] = src.astype(acqs[i].data.dtype, copy=False)
        acqs[i].data[:] = tmp

        trim_lengths[i] = n

    # step 2: optionally resize acquisitions
    if mode == "trim":
        for i in acq_idx:
            n_new = trim_lengths[i]
            nchan = acqs[i].data.shape[0]

            data_trim = acqs[i].data[:, :n_new].copy()

            has_traj = hasattr(acqs[i], "traj") and (acqs[i].traj is not None)
            if has_traj:
                traj_trim = acqs[i].traj[:n_new, :].copy()
                traj_dim = traj_trim.shape[1]
            else:
                traj_trim = None
                traj_dim = 0

            acqs[i].resize(
                number_of_samples=n_new,
                active_channels=nchan,
                trajectory_dimensions=traj_dim
            )

            acqs[i].data[:] = data_trim
            if traj_trim is not None:
                acqs[i].traj[:] = traj_trim

def add_bonus_idx_to_user_params(userParams, bonus_idx):
    """
    Store bonus_idx in MRD header user parameters.

    Saves:
      - bonus_spectra_present : 0/1
      - bonus_spectra_count   : integer count
      - bonus_spectra_indices : comma-separated string of indices
    """
    if bonus_idx is None:
        bonus_idx = []

    bonus_idx = [int(i) for i in bonus_idx]

    bonus_present = mrd.xsd.userParameterLongType('bonus_spectra_present')
    bonus_present.value = int(len(bonus_idx) > 0)
    userParams.userParameterLong.insert(0, bonus_present)

    bonus_count = mrd.xsd.userParameterLongType('bonus_spectra_count')
    bonus_count.value = len(bonus_idx)
    userParams.userParameterLong.insert(0, bonus_count)

    bonus_indices = mrd.xsd.userParameterStringType('bonus_spectra_indices')
    bonus_indices.value = ",".join(str(i) for i in bonus_idx)
    userParams.userParameterString.insert(0, bonus_indices)

    return userParams

def reorder_crds_to_scanner_labels(acqs, crds, rep_to_use=None, set_to_use=None):
    """
    Convert crds from acquisition-order profile columns to scanner ky-label columns.

    Input:
        crds shape = [n_interleaves, n_profiles, n_samples, 3]
        acqs = list of acquisitions in original acquisition order

    Output:
        crds_out shape = same as crds, but column index now matches acq.idx.kspace_encode_step_1
    """
    n_interleaves, n_profiles = crds.shape[0], crds.shape[1]

    # Build ky appearance order from one repetition / one echo only
    ky_order = []
    seen = set()

    for a in acqs:
        if a is None:
            continue
        if rep_to_use is not None and int(a.idx.repetition) != int(rep_to_use):
            continue
        if set_to_use is not None and int(a.idx.set) != int(set_to_use):
            continue

        ky = int(a.idx.kspace_encode_step_1)
        kz = int(a.idx.kspace_encode_step_2)

        # only take first interleave block entry for each ky
        if kz == 0 and ky not in seen:
            ky_order.append(ky)
            seen.add(ky)

    # fallback if kz indexing is not zero-based in this reader
    if len(ky_order) != n_profiles:
        ky_order = []
        seen = set()
        min_kz = min(int(a.idx.kspace_encode_step_2) for a in acqs if a is not None)
        for a in acqs:
            if a is None:
                continue
            if rep_to_use is not None and int(a.idx.repetition) != int(rep_to_use):
                continue
            if set_to_use is not None and int(a.idx.set) != int(set_to_use):
                continue

            ky = int(a.idx.kspace_encode_step_1)
            kz = int(a.idx.kspace_encode_step_2)

            if kz == min_kz and ky not in seen:
                ky_order.append(ky)
                seen.add(ky)

    if len(ky_order) != n_profiles:
        raise ValueError(
            f"Could not build ky permutation correctly: got {len(ky_order)} profiles, expected {n_profiles}"
        )

    ky_order = np.asarray(ky_order, dtype=int)

    # handle possible 1-based labels
    if ky_order.min() == 1 and ky_order.max() == n_profiles:
        ky_order = ky_order - 1

    if ky_order.min() != 0 or ky_order.max() != n_profiles - 1:
        raise ValueError(
            f"Unexpected ky label range: min={ky_order.min()}, max={ky_order.max()}, expected 0..{n_profiles-1}"
        )

    crds_out = np.empty_like(crds)

    # column j in crds = j-th acquired ky block
    # move it into the scanner's ky label column
    for j, ky_label in enumerate(ky_order):
        crds_out[:, ky_label, :, :] = crds[:, j, :, :]

    return crds_out, ky_order

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
    n_acq_use = n_acq_total
    bonus_idx = []
    bonus_info = {}
    if data_set_config.data_type == DataType.DIXON:     
        #if data_set_config.gas_contam_removal or data_set_config.exclude_bonus_spec:
        acqs_init      = [None] * n_acq_total
        for acqnum in range(n_acq_total):
            acq_temp = dset.read_acquisition(acqnum)
            acqs_init[acqnum] = acq_temp
        bonus_idx, bonus_info = find_odd_spectra_indices(acqs_init)
        print('bonus_idx = ', bonus_idx)
    
        n_acq_use = n_acq_total - (len(bonus_idx) if data_set_config.exclude_bonus_spec else 0)
        print("n_acq_total =", n_acq_total, "n_acq_use =", n_acq_use,
            "exclude_bonus_spec =", data_set_config.exclude_bonus_spec)

    if data_set_config.data_type == DataType.DIXON: 
        # If requested, create a new MRD file that excludes the bonus_spec acquisition
        if data_set_config.exclude_bonus_spec and not data_set_config.gas_contam_removal:

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
            
    debug_mode = data_set_config.debug_mode 
    
    # Get path to sin file trajectories if necessary
    if data_set_config.data_type == DataType.CALIBRATION: 
        data_set_config.ext_traj == False
    if data_set_config.ext_traj == True:
        if traj_file == None:
            # Directory containing the running script
            script_dir = Path(__file__).resolve().parent
            repo_root = script_dir.parent

            if 'FLORET'.lower() in rls.header['sin']['scan_name'][0][0].lower():
                ext_coords = repo_root / "resources" / "Polarean_Xenon_FLORET_Dixon_20251205-NoSpectra.sin"
            else:
                ext_coords = repo_root / "resources" / "20221025_144528_DukeIPF_Gas_Exchange.sin"
            
            if data_set_config.institution == 'CCHMC':
                ext_coords = repo_root / "resources" / "CCHMC_Dissolved_Xe_20191008 - 3T-T1.sin"

            if data_set_config.institution == 'Polarean':
                ext_coords = repo_root / "resources" / "Polarean_Xenon_Radial_Dixon_20260211_NoSpec.sin"
                
            # Safety check
            if not ext_coords.exists():
                ext_coords = filedialog.askopenfilename(title='Select trajectory .sin file', filetypes=[
                    ("Trajectory .sin file", "*.sin")], initialdir=outDir)
        else:
            ext_coords = traj_file
            
        inputDataTraj = rp.PhilipsData(ext_coords)
        inputDataTraj.trajtype = data_set_config.trajorder
        print('config_traj_type = ', inputDataTraj.trajtype)
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
            
        # Haltoned Spiral, multiple interleaves: remap columns from acquisition order -> scanner ky labels
        if data_set_config.data_type == DataType.DIXON: 
            n_interleaves = crds.shape[0]
            if data_set_config.trajorder == 2 and n_interleaves > 1:
                acqs_for_perm = [dset.read_acquisition(i) for i in range(dset.number_of_acquisitions())]

                # use first repetition only to avoid duplicate gas/diss blocks
                rep0 = min(int(a.idx.repetition) for a in acqs_for_perm if a is not None)
                crds, ky_order = reorder_crds_to_scanner_labels(acqs_for_perm, crds, rep_to_use=rep0)
                if debug_mode:
                    print("ky acquisition order -> ky label:")
                    print(ky_order[:50])

                if 'crds_flyback' in locals():
                    crds_flyback, _ = reorder_crds_to_scanner_labels(acqs_for_perm, crds_flyback, rep_to_use=rep0)            
            
    if debug_mode:      
        # save coords as .mat file
        save_dir = os.path.dirname(data_file)
        os.makedirs(save_dir, exist_ok=True)
        save_path = os.path.join(save_dir, "crds.mat")
        # save
        savemat(save_path, {'crds': crds.astype('float64')})
        print(f"Saved to: {save_path}")                    
  
    if data_set_config.data_type == DataType.DIXON and debug_mode:        
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
    
    # Save bonus spectra indices in header
    if data_set_config.data_type == DataType.DIXON:
        if data_set_config.gas_contam_removal and not data_set_config.exclude_bonus_spec:    
            userParams = add_bonus_idx_to_user_params(userParams, bonus_idx)
    
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
    if data_set_config.data_type in (DataType.DIXON, DataType.UTE):
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


    # Read/store all data/labels for ease/simplicity
    # ---------------- Gas contamination removal ----------------
        
    if data_set_config.gas_contam_removal:
        # -----------------------------
        # Handle bonus spectra
        # -----------------------------
        if bonus_idx is None or len(bonus_idx) == 0:
            print("No bonus spectra found by size; skipping gas contamination removal.")
            bonus_idx = []
            bonus_acqs = []
        else:
            bonus_idx = [int(i) for i in bonus_idx]
            bonus_acqs = [acqs_init[i] for i in bonus_idx]

        bonus_set = set(bonus_idx)

        # Keep only non-bonus acquisitions if requested
        if data_set_config.exclude_bonus_spec:
            keep_idx = [i for i in range(n_acq_total) if i not in bonus_set]
        else:
            keep_idx = list(range(n_acq_total))

        n_acq_use = len(keep_idx)

        print(
            "n_acq_total =", n_acq_total,
            "n_acq_use =", n_acq_use,
            "exclude_bonus_spec =", data_set_config.exclude_bonus_spec
        )

        if n_acq_use == 0:
            raise ValueError("No acquisitions remain after excluding bonus spectra.")

        # -----------------------------
        # Find oversampling factor
        # -----------------------------
        acq_temp = dset.read_acquisition(keep_idx[3])
        data0 = acq_temp.data
        Npoints = data0.shape[-1]
        traj_Npoints = crds.shape[2]
        OvsFactor = int(Npoints / traj_Npoints)

        # -----------------------------
        # Read kept acquisitions only
        # -----------------------------
        acqs = []
        data = []
        contrasts = []
        sets = []
        kzs = []
        kys = []
        kxs = []

        for acqnum in keep_idx:
            acq_temp = dset.read_acquisition(acqnum)
            acqs.append(acq_temp)

            if OvsFactor > 1:
                data02 = movmean(acq_temp.data, OvsFactor, axis=1)
                data02 = data02[:, ::OvsFactor]
                data.append(data02)
                kxs.append(data02.shape[1])
            else:
                data.append(acq_temp.data)
                kxs.append(acq_temp.data.shape[1])

            contrasts.append(acq_temp.idx.contrast)
            sets.append(acq_temp.idx.set)
            kzs.append(acq_temp.idx.kspace_encode_step_2)
            kys.append(acq_temp.idx.kspace_encode_step_1)

        # -----------------------------
        # Bonus FIDs
        # -----------------------------
        if len(bonus_acqs) > 0:
            bonus_fids = [np.squeeze(acq.data) for acq in bonus_acqs]
            bonus_fid = np.array(bonus_fids)

            if debug_mode:
                plot_complex_spectra_panels(bonus_fid, do_fft=False, suptitle='Bonus Spectra')

            print('bonus_fid size =', bonus_fid[0].size)
            dwell_s = float(dwell.value) * 1e-6

            n_fids = len(bonus_fid)
            inst = data_set_config.institution

            if n_fids == 6: # and inst == 'CCHMC':
                PreDissolvedFID  = bonus_fid[0]
                PostDissolvedFID = bonus_fid[3]
                PreGasFID        = bonus_fid[1]
                PostGasFID       = bonus_fid[4]

            elif n_fids == 2: # and inst == 'Polarean':
                zeros_fid = np.zeros_like(bonus_fid[0])
                PreDissolvedFID  = zeros_fid
                PostDissolvedFID = bonus_fid[0]
                PreGasFID        = zeros_fid
                PostGasFID       = bonus_fid[1]

            else:
                raise ValueError(
                    f"Unsupported combination: n_fids={n_fids}, institution={inst}"
                )
        else:
            bonus_fid = None
            PreDissolvedFID = None
            PostDissolvedFID = None
            PreGasFID = None
            PostGasFID = None

        # -----------------------------
        # K-space split
        # -----------------------------
        gas_rep  = data_set_config.contrast_order.index(1)
        diss_rep = data_set_config.contrast_order.index(2)

        # acqs already contains only the acquisitions you decided to keep,
        # so do NOT filter with bonus_set again here
        gas_acq_idx = [i for i, a in enumerate(acqs) if a is not None and a.idx.repetition == gas_rep]
        diss_acq_idx = [i for i, a in enumerate(acqs) if a is not None and a.idx.repetition == diss_rep]       
        gas_acq_idx = gas_acq_idx[1:-1]
        diss_acq_idx = diss_acq_idx[1:-1]

        if len(gas_acq_idx) == 0:
            raise ValueError("No gas acquisitions found after filtering.")

        if len(diss_acq_idx) == 0:
            raise ValueError("No dissolved acquisitions found after filtering.")

        GasKSpaceInit = np.stack(
            [np.squeeze(acqs[i].data).reshape(-1) for i in gas_acq_idx],
            axis=1
        )

        DissolvedKSpaceInit = np.stack(
            [np.squeeze(acqs[i].data).reshape(-1) for i in diss_acq_idx],
            axis=1
        )

        #OvsFactor = 1
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

        if debug_mode:
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
        # --- 1) Fit prepended dissolved (Pre) ---
        PrependedDissolvedNMRFit = NMR_TimeFit(
            ydata=PreDissolvedFID,
            tdata=time,
            area=area_guess,
            freq=freq_guess,
            fwhmL=fwhmL_guess,
            fwhmG=fwhmG_guess,
            phase=phase_guess,
            line_broadening=0,
            zeropad_size=time.size,
            method="voigt",
        )
        PrependedDissolvedNMRFit = make_timefit_x0_feasible(PrependedDissolvedNMRFit, lb, ub)
        PrependedDissolvedNMRFit.fit_time_signal_residual(bounds=(lb, ub))
     
        # -------------------- 2) Fit appended dissolved (Post), seeded by Pre-fit --------------------
        if PreDissolvedFID is not None and np.any(np.abs(PreDissolvedFID) > 0):
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
        else:
            AppendedDissolvedNMRFit = NMR_TimeFit(
                ydata=PostDissolvedFID,
                tdata=time,
                area=area_guess,
                freq=freq_guess,
                fwhmL=fwhmL_guess,
                fwhmG=fwhmG_guess,
                phase=phase_guess,
                line_broadening=0,
                zeropad_size=time.size,
                method="voigt",
            )            
            
        # Optional: initial fit (unbounded) like MATLAB
        #AppendedDissolvedNMRFit.fit_time_signal_residual(bounds=(-np.inf, np.inf))

        # Critical: make x0 feasible for bounded refit (SciPy requirement)
        AppendedDissolvedNMRFit = make_timefit_x0_feasible(AppendedDissolvedNMRFit, lb, ub)

        # Bounded refit (MATLAB "setBounds + refit")
        AppendedDissolvedNMRFit.fit_time_signal_residual(bounds=(lb, ub))
        
        # -------------------- 3)  dwell_s * fftshift(fft(calcComponentTimeDomainSignal(time),[],1),1) --------------------
        comp_td = AppendedDissolvedNMRFit.calcComponentTimeDomainSignal(time)  # (Nt, 3)

        # FFT each component along time axis (axis=0), then fftshift along same axis
        fit_components = dwell_s * np.fft.fftshift(np.fft.fft(comp_td, axis=0), axes=0)  # (Nf, 3)

        # Sum of components (total fit spectrum)
        fit_sum = np.sum(fit_components, axis=1)  # (Nf,)

        phase = AppendedDissolvedNMRFit.phase[2]
        area  = AppendedDissolvedNMRFit.area[2]

        if debug_mode:
            f = AppendedDissolvedNMRFit.f
            spec = AppendedDissolvedNMRFit.spectral_signal  # measured spectrum

            plt.figure(figsize=(16, 9), facecolor="white")

            # ---------------- Data ----------------
            plt.plot(f, np.abs(spec), 'b', label='Spectrum - Magnitude')
            plt.plot(f, np.real(spec), 'r', label='Spectrum - Real')
            plt.plot(f, np.imag(spec), color=(0, 0.6, 0.2), label='Spectrum - Imaginary')

            # ---------------- Fits (sum) ----------------
            plt.plot(f, np.abs(fit_sum), color=(0, 0, 1, 0.33), linewidth=3, label='Fit - Magnitude')
            plt.plot(f, np.real(fit_sum), color=(1, 0, 0, 0.33), linewidth=3, label='Fit - Real')
            plt.plot(f, np.imag(fit_sum), color=(0, 0.6, 0.2, 0.33), linewidth=3, label='Fit - Imaginary')

            # ---------------- Components (Real only) ----------------
            plt.fill_between(f, 0, np.real(fit_components[:, 0]), color='r', alpha=0.33, label='RBC - Real')
            plt.fill_between(f, 0, np.real(fit_components[:, 1]), color='b', alpha=0.33, label='Barrier - Real')
            plt.fill_between(f, 0, np.real(fit_components[:, 2]), color=(0, 0.6, 0.2), alpha=0.33, label='Gas - Real')

            plt.legend(loc='best', ncol=3, frameon=False)
            plt.xlim([-8000, 2000])
            plt.gca().invert_xaxis()
            plt.xticks(fontsize=22)
            plt.yticks(fontsize=22)
            
            leg = plt.legend(loc='best', ncol=3, frameon=False)
            for text in leg.get_texts():
                text.set_fontsize(18)

            ratio = AppendedDissolvedNMRFit.area[0] / AppendedDissolvedNMRFit.area[1]
            plt.title(f'Dissolved Phase Spectra and Fit: RBC/Barrier = {ratio:.3f}', fontsize=22)
            plt.xlabel('Frequency (Hz)', fontsize=20)
            plt.ylabel('NMR Signal (a.u.)', fontsize=20)

            plt.tight_layout()
            plt.show()

            print('Fitting Spectrum Completed.')

        freq_jump = data_set_config.xe_dissolved_offset_ppm * centFreq.value / 1e6

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

        if debug_mode:
            # Equivalent to size(...,2)
            XePulses = DissolvedKSpaceInit.shape[1]
            SS_ind = 60
            # Create figure
            fig = plt.figure(figsize=(16, 9))
            fig.patch.set_facecolor('white')

            # Subplot
            ax = fig.add_subplot(1, 2, 1)
            ax.tick_params(labelsize=18)

            # --- Data ---
            # Dissolved (original)
            ax.plot(np.abs(DissolvedKSpaceInit[0, :]),
                    color=(0.5, 0, 0), linewidth=1, label='Dissolved')

            # Corrected dissolved (MATLAB used XePulses as x, but that plots a single point;
            # assuming you meant full trace)
            ax.plot(np.abs(DissKSpace_corr[0, :]),
                    color=(1, 0, 0), linewidth=1, label='Corrected Dissolved')

            # Gas
            ax.plot(np.abs(GasKSpaceInit[0, :]),
                    color=(0, 0, 1), linewidth=1, label='Gas')

            # --- Steady-state shading ---
            yl = ax.get_ylim()
            xl = ax.get_xlim()

            ax.fill([0, SS_ind, SS_ind, 0],
                    [yl[1], yl[1], yl[0], yl[0]],
                    color=(0.66, 0.66, 0.66),
                    alpha=1.0,
                    label='Steady State Discard',
                    zorder=0)

            # Restore limits (like MATLAB)
            ax.set_xlim(xl)
            ax.set_ylim(yl)

            # --- Labels ---
            ax.set_title('Xe Signal Intensity Dynamics')
            ax.set_xlabel('Projection')
            ax.set_ylabel('Signal Intensity')
            ax.legend(loc='best', ncol=2, fontsize=12)
            ax.grid(False)

            plt.tight_layout()
            plt.show()    
            
        # e) write corrected data back into dissolved acquisitions
        mode = "trim" # mode: 'pad' (do nothing) or 'trim'
        # dissolved
        write_kspace_to_acqs(
            acqs=acqs,
            acq_idx=diss_acq_idx,
            kspace_mat=DissKSpace_corr,
            mode=mode
        )
        # gas
        write_kspace_to_acqs(
            acqs=acqs,
            acq_idx=gas_acq_idx,
            kspace_mat=GasKSpaceInit,
            mode=mode
        )        
        print("Gas contamination removal applied.")

    else:
        acqs      = [None] * n_acq_use
        data      = [None] * n_acq_use
        contrasts = [None] * n_acq_use
        sets      = [None] * n_acq_use
        kzs       = [None] * n_acq_use
        kys       = [None] * n_acq_use
        kxs       = [None] * n_acq_use

        if data_set_config.data_type == DataType.DIXON:
            acq_temp = dset.read_acquisition(10)
            data0 = acq_temp.data
            Npoints = data0.shape[-1]
            traj_Npoints = crds.shape[2]
            OvsFactor = int(Npoints / traj_Npoints)

            for acqnum in range(len(acqs)):
                acq_temp = dset.read_acquisition(acqnum)
                data_i = acq_temp.data

                if OvsFactor > 1:
                    data_ds = movmean(data_i, OvsFactor, axis=1)
                    data_ds = data_ds[:, ::OvsFactor]

                    nchan = data_ds.shape[0]
                    nsamp_new = data_ds.shape[1]

                    has_traj = hasattr(acq_temp, "traj") and (acq_temp.traj is not None)
                    if has_traj:
                        traj_old = acq_temp.traj.copy()
                        traj_dim = traj_old.shape[1]
                    else:
                        traj_old = None
                        traj_dim = 0

                    acq_temp.resize(
                        number_of_samples=nsamp_new,
                        active_channels=nchan,
                        trajectory_dimensions=traj_dim
                    )

                    acq_temp.data[:] = data_ds.astype(acq_temp.data.dtype, copy=False)

                    if traj_old is not None:
                        traj_trim = traj_old[:nsamp_new, :]
                        acq_temp.traj[:] = traj_trim.astype(acq_temp.traj.dtype, copy=False)

                acqs[acqnum] = acq_temp
                data[acqnum] = acq_temp.data.copy()
                contrasts[acqnum] = acq_temp.idx.contrast
                sets[acqnum] = acq_temp.idx.set
                kzs[acqnum] = acq_temp.idx.kspace_encode_step_2
                kys[acqnum] = acq_temp.idx.kspace_encode_step_1
                kxs[acqnum] = acq_temp.number_of_samples
        else:
            for acqnum in range(n_acq_use):
                acq_temp = dset.read_acquisition(acqnum)
                acqs[acqnum] = acq_temp
                data[acqnum] = acq_temp.data.copy()
                contrasts[acqnum] = acq_temp.idx.contrast
                sets[acqnum] = acq_temp.idx.set
                kzs[acqnum] = acq_temp.idx.kspace_encode_step_2
                kys[acqnum] = acq_temp.idx.kspace_encode_step_1
                kxs[acqnum] = acq_temp.number_of_samples    
    # -----------------------------------------------------------

    # update data and headers
    # Version   |   
    # Philips:  | location | average | extr1   | mix      | card  | dynamic    | echo     | kz | ky | kx
    # MRD:      | slice    | average | segment | set      | phase | repetition | contrast | kz | ky | kx
    # XeCTCMRD: | N/A      | N/A     | N/A     | spec/img | N/A   | contrast   | set      | kz | ky | kx
    
    for acqnum in range(len(acqs)):
        acq_temp = acqs[acqnum]
        if debug_mode:
            if acqnum < 20:
                print(
                    f"acq {acqnum}: "
                    f"contrast={acq_temp.idx.contrast}, "
                    f"rep={acq_temp.idx.repetition}, "
                    f"kz={acq_temp.idx.kspace_encode_step_2}, "
                    f"ky={acq_temp.idx.kspace_encode_step_1}"
                )
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
                if debug_mode:
                    if acqnum < 5:
                        print("traj first point:", traj[0, :])
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
        
    if data_set_config.data_type == DataType.DIXON:
        if data_set_config.gas_contam_removal and data_set_config.exclude_bonus_spec:
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
                raw_file=args.raw_file, 
                traj_file=args.traj_file)


