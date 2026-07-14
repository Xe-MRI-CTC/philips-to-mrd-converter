def read_ismrmrd_file(filepath, lean=False):
    import os
    import ismrmrd
    import ismrmrd.xsd
    import numpy as np

    if not os.path.isfile(filepath):
        print("%s is not a valid file" % filepath)
        raise SystemExit
    dset = ismrmrd.Dataset(filepath, "dataset", create_if_needed=False)

    header = ismrmrd.xsd.CreateFromDocument(dset.read_xml_header())
    enc = header.encoding[0]

    # Matrix size
    eNx = enc.encodedSpace.matrixSize.x
    eNy = enc.encodedSpace.matrixSize.y
    eNz = enc.encodedSpace.matrixSize.z
    rNx = enc.reconSpace.matrixSize.x
    rNy = enc.reconSpace.matrixSize.y
    rNz = enc.reconSpace.matrixSize.z

    # Field of View
    eFOVx = enc.encodedSpace.fieldOfView_mm.x
    eFOVy = enc.encodedSpace.fieldOfView_mm.y
    eFOVz = enc.encodedSpace.fieldOfView_mm.z
    rFOVx = enc.reconSpace.fieldOfView_mm.x
    rFOVy = enc.reconSpace.fieldOfView_mm.y
    rFOVz = enc.reconSpace.fieldOfView_mm.z

    # Number of Slices, Reps, Contrasts, etc.
    ncoils = header.acquisitionSystemInformation.receiverChannels
    if enc.encodingLimits.slice is not None:
        nslices = enc.encodingLimits.slice.maximum + 1
    else:
        nslices = 1

    if enc.encodingLimits.repetition is not None:
        nreps = enc.encodingLimits.repetition.maximum + 1
    else:
        nreps = 1

    if enc.encodingLimits.contrast is not None:
        ncontrasts = enc.encodingLimits.contrast.maximum + 1
    else:
        ncontrasts = 1

    # TODO loop through the acquisitions looking for noise scans
    firstacq = 0
    for acqnum in range(dset.number_of_acquisitions()):
        acq = dset.read_acquisition(acqnum)

        # TODO: Currently ignoring noise scans
        if acq.isFlagSet(ismrmrd.ACQ_IS_NOISE_MEASUREMENT):
            print("Found noise scan at acq ", acqnum)
            continue
        else:
            firstacq = acqnum
            print("Imaging acquisition starts acq ", acqnum)
            break

    # Initialiaze a storage array
    # Dimensions
    nDims = np.shape(acq.traj)[1]
    if lean:
        all_data = np.zeros((1, 1, nslices, ncoils, eNz, eNy, eNx), dtype=np.complex64)
        all_traj = np.zeros((1, 1, nslices, eNz, eNy, eNx, nDims), dtype=np.float64)
    else:
        all_data = np.zeros((nreps, ncontrasts, nslices, ncoils, eNz, eNy, eNx), dtype=np.complex64)
        all_traj = np.zeros((nreps, ncontrasts, nslices, eNz, eNy, eNx, nDims), dtype=np.float64)

    # Loop through the rest of the acquisitions and stuff
    for acqnum in range(firstacq, dset.number_of_acquisitions()):
        acq = dset.read_acquisition(acqnum)

        # Stuff into the buffer
        rep = acq.idx.repetition
        contrast = acq.idx.contrast
        slice = acq.idx.slice
        y = acq.idx.kspace_encode_step_1
        z = acq.idx.kspace_encode_step_2

        if lean:
            if rep > enc.encodingLimits.repetition.minimum or contrast > enc.encodingLimits.contrast.minimum:
                continue

        all_data[rep, contrast, slice, :, z, y, :] = acq.data
        all_traj[rep, contrast, slice, z, y, :, :] = acq.traj

    return {"data": all_data, "traj": all_traj, "header": header}


def view_image_montage(image, axis="z", rows=None, cols=None, cmap="gray", title=None):
    """
    Create a montage of 2D slices from a 3D volume.

    Parameters:
    volume: numpy.ndarray
        3D volume array
    axis: str, optional
        Axis along which to slice ('x', 'y', 'z')
    rows, cols: int, optional
        Number of rows and columns in the montage grid. If not provided,
        will attempt to create a square-ish grid
    cmap: str, optional
        Colormap for matplotlib
    title: str, optional
        Title for the plot
    """
    import numpy as np
    import matplotlib.pyplot as plt

    # Map axis to axis index
    axis_map = {"x": 0, "y": 1, "z": 2}
    if axis not in axis_map:
        raise ValueError("Axis must be 'x', 'y', or 'z'")

    axis_idx = axis_map[axis]

    # Get the number of slices along the chosen axis
    n_slices = image.shape[axis_idx]

    # Determine grid dimensions if not provided
    if rows is None and cols is None:
        cols = int(np.ceil(np.sqrt(n_slices)))
        rows = int(np.ceil(n_slices / cols))
    elif rows is None:
        rows = int(np.ceil(n_slices / cols))
    elif cols is None:
        cols = int(np.ceil(n_slices / rows))

    # Create the montage figure
    fig, axes = plt.subplots(rows, cols, figsize=(cols * 2, rows * 2))

    # Handle case where only one subplot
    if rows == 1 and cols == 1:
        axes = [axes]
    else:
        axes = axes.flatten()

    # Loop through slices and plot
    for i, ax in enumerate(axes):
        if i < n_slices:
            # Extract slice based on axis
            if axis_idx == 0:
                slice_img = image[i, :, :]
            elif axis_idx == 1:
                slice_img = image[:, i, :]
            else:  # axis_idx == 2
                slice_img = image[:, :, i]

            ax.imshow(slice_img, cmap=cmap)
            ax.set_title(f"Slice {i+1}/{n_slices}")
            ax.axis("off")
        else:
            ax.axis("off")

    if title:
        fig.suptitle(title)

    plt.tight_layout()
    plt.show()


def recon_noncart_mrd(mrd_dict):
    import sigpy as sp
    import sigpy.mri as mr
    import numpy as np

    traj = mrd_dict["traj"]  # allocate
    n_dims = np.size(traj, -1)

    if n_dims == 3:
        recon_mtx = [
            mrd_dict["header"].encoding[0].reconSpace.matrixSize.x,
            mrd_dict["header"].encoding[0].reconSpace.matrixSize.y,
            mrd_dict["header"].encoding[0].reconSpace.matrixSize.z,
        ]
        traj[..., 0] = mrd_dict["traj"][..., 0] * recon_mtx[0]
        traj[..., 1] = mrd_dict["traj"][..., 1] * recon_mtx[1]
        traj[..., 2] = mrd_dict["traj"][..., 2] * recon_mtx[2]
        k_data = np.squeeze(mrd_dict["data"])
    else:
        recon_mtx = [
            mrd_dict["header"].encoding[0].reconSpace.matrixSize.x,
            mrd_dict["header"].encoding[0].reconSpace.matrixSize.y,
        ]
        n_slices = mrd_dict["header"].encoding[0].encodingLimits.slice.maximum + 1
        n_projs = mrd_dict["header"].encoding[0].encodedSpace.matrixSize.y
        traj[..., 0] = mrd_dict["traj"][..., 0] * recon_mtx[0]
        traj[..., 1] = mrd_dict["traj"][..., 1] * recon_mtx[1]
        traj = traj[0, 0, 0, 0, :, :, :]
        k_data = np.squeeze(mrd_dict["data"])

    F = sp.linop.NUFFT(recon_mtx, coord=traj, oversamp=1.25, width=4, toeplitz=True)

    dcf = mr.pipe_menon_dcf(traj, img_shape=recon_mtx, beta=8, width=3)
    dcf /= np.max(dcf)
    D = sp.linop.Multiply(F.oshape, dcf**0.5)

    F_dcf = D * F

    if n_dims == 3:
        k_data_dcf = k_data * dcf**0.5
        image_nufft = abs(F_dcf.H * k_data_dcf)
    else:
        image_nufft = np.zeros([n_slices, recon_mtx[0], recon_mtx[1]])
        for s in range(0, n_slices):
            k_data_temp = k_data[s, :, :]
            k_data_dcf = k_data_temp * dcf**0.5
            image_nufft[s, :, :] = abs(F_dcf.H * k_data_dcf)

    return image_nufft


def recon_cart_mrd(mrd_dict):
    import numpy as np

    dims = 2 if mrd_dict["header"].encoding[0].reconSpace.matrixSize.z == 1 else 3
    k_data = mrd_dict["data"]
    n_reps, n_contr, n_slices, n_channels, nZ, nY, nX = np.shape(k_data)

    if n_reps > 1 or n_contr > 1:
        print("Only reconsrtucting first rep and echo")

    k_data = k_data[0, 0, ::]  # remove rep and contr dims
    k_data = np.transpose(k_data, (1, 0, 2, 3, 4))  # move channels to outer dim
    k_data = np.squeeze(k_data)  # either slice or nz should be 1 and removed

    if dims == 2:
        axes = (-2, -1)
    else:
        axes = (-3, -2, -1)

    k_data = np.fft.ifftshift(k_data, axes=axes)
    image = np.fft.ifftn(k_data, axes=axes)
    image = np.fft.fftshift(image, axes=(-1))
    image = np.abs(image)

    if np.shape(k_data)[0] > 1 and k_data.ndim == 4:  # if multiple channels
        image = np.sqrt(np.sum(image**2, axis=0))
    elif k_data.ndim == 4:
        image = image[0, ::]

    x_chop = int((nX - mrd_dict["header"].encoding[0].reconSpace.matrixSize.x) / 2)
    y_chop = int((nY - mrd_dict["header"].encoding[0].reconSpace.matrixSize.y) / 2)
    image = image[:, y_chop:-y_chop, x_chop:-x_chop]

    return image


def recon_noncart_rp(rp_results):
    import sigpy as sp
    import sigpy.mri as mr
    import numpy as np

    traj = rp_results.coords  # allocate
    n_dims = int(rp_results.header['sin']['encoding_dimensions'][0][0])

    k_data = rp_results.data
    if rp_results.header['headerType'] == 'lab-sin':
        # FROM PhilipsData:     data_string = np.array(['chan', 'mix', 'dyn', 'card', 'echo', 'row',
        #                                               'extra', 'loc', 'e3', 'meas', 'e2', 'e1', 'samp'])
        k_data = k_data[:,0,0,0,0,0,0,:,0,0,:,:,:] #loc, kz, ky, samp
        k_data = np.squeeze(k_data)  # either slice or nz should be 1 and removed

    if n_dims == 3:
        recon_mtx = [
            int(rp_results.header['sin']['recon_resolutions'][0][0]),
            int(rp_results.header['sin']['recon_resolutions'][0][1]),
            int(rp_results.header['sin']['recon_resolutions'][0][2]),
        ]
        traj[..., 0] = rp_results.coords[..., 0] * recon_mtx[0]
        traj[..., 1] = rp_results.coords[..., 1] * recon_mtx[1]
        traj[..., 2] = rp_results.coords[..., 2] * recon_mtx[2]
    else:
        recon_mtx = [
            int(rp_results.header['sin']['recon_resolutions'][0][0]),
            int(rp_results.header['sin']['recon_resolutions'][0][1]),
        ]
        n_slices = int(rp_results.header['sin']['nr_locations'][0][0])
        n_projs = int(rp_results.header['sin']['non_cart_max_encoding_nrs'][0][1]) + 1
        traj[..., 0] = rp_results.coords[..., 0] * recon_mtx[0]
        traj[..., 1] = rp_results.coords[..., 1] * recon_mtx[1]
        traj = np.delete(traj, 2, axis=-1)

    F = sp.linop.NUFFT(recon_mtx, coord=traj, oversamp=1.25, width=4, toeplitz=True)

    dcf = mr.pipe_menon_dcf(traj, img_shape=recon_mtx, beta=8, width=3)
    dcf /= np.max(dcf)
    D = sp.linop.Multiply(F.oshape, dcf**0.5)

    F_dcf = D * F

    if n_dims == 3:
        k_data_dcf = k_data * dcf**0.5
        image_nufft = abs(F_dcf.H * k_data_dcf)
    else:
        image_nufft = np.zeros([n_slices, recon_mtx[0], recon_mtx[1]])
        for s in range(0, n_slices):
            k_data_temp = k_data[s, :, :]
            k_data_dcf = k_data_temp * dcf**0.5
            image_nufft[s, :, :] = abs(F_dcf.H * k_data_dcf)

    return image_nufft


def recon_cart_rp(rp_results):
    import numpy as np

    k_data = rp_results.data
    if rp_results.header['headerType'] == 'lab-sin':
        dims = int(rp_results.header['sin']['encoding_dimensions'][0][0])
        # FROM PhilipsData:     data_string = np.array(['chan', 'mix', 'dyn', 'card', 'echo', 'row',
        #                                               'extra', 'loc', 'e3', 'meas', 'e2', 'e1', 'samp'])
        k_data = k_data[:,0,0,0,0,0,0,:,0,0,:,:,:] #loc, kz, ky, samp
        k_data = np.squeeze(k_data)  # either slice or nz should be 1 and removed

    if dims == 2:
        axes = (-2, -1)
    else:
        axes = (-3, -2, -1)

    k_data = np.fft.ifftshift(k_data, axes=axes)
    image = np.fft.ifftn(k_data, axes=axes)
    image = np.fft.fftshift(image, axes=(-1))
    image = np.abs(image)

    if np.shape(k_data)[0] > 1 and k_data.ndim == 4:  # if multiple channels
        image = np.sqrt(np.sum(image**2, axis=0))
    elif k_data.ndim == 4:
        image = image[0, ::]
    
    nX = image.shape[-1]
    nY = image.shape[-2]
    x_chop = int((nX - int(rp_results.header['sin']['scan_resolutions'][0][0])) / 2)
    y_chop = int((nY - int(rp_results.header['sin']['scan_resolutions'][0][1])) / 2)
    image = image[:, y_chop:-y_chop, x_chop:-x_chop]

    return image