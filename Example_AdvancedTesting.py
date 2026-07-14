from Scripts import XeGasExchange2XeCTCMRD
import readphilips.ReadPhilips as rp
import numpy as np
from philips2mrd.functions import mrd_compare as mrc
from philips2mrd.functions import mrd_recon as mrr
import philips2mrd.philips2mrd as p2m
from tests import datasets
import matplotlib.pyplot as plt

# -----------SELECT DATA-----------#
test_data = datasets.getTestData()
data_number = 1

data_type = test_data[data_number]["type"]
data_name = test_data[data_number]["name"]
data_file = test_data[data_number]["dl"]
raw_file = test_data[data_number]["rls"]
traj_file = test_data[data_number]["traj"]


# -----------BASIC READ PHILIPS-----------#
result = rp.PhilipsData(raw_file)
result.compute()


# -----------READ PHILIPS RECONSTRUCT-----------#
traj_type = 0
try:
    traj_type = int(result.header["sin"]["k_space_traj_type"][0][0])
except Exception:
    traj_type = 0

# Basic data views
if traj_type == 0:
    readouts = np.squeeze(result.data)
    readouts = readouts.reshape(-1, readouts.shape[-1])
    fig, ax = plt.subplots()
    ax.plot(range(np.shape(readouts)[-1]), np.transpose(readouts, axes=(1,0)))
    plt.show(block=False)
else:
    k0 = result.data[..., 0]
    k0 = abs(k0.reshape(-1))
    fig, ax = plt.subplots()
    ax.plot(range(len(k0)), k0, "-o")
    plt.show(block=False)

    kx_max = result.coords[..., -1, 0]
    kx_max = abs(kx_max.reshape(-1))
    fig, ax = plt.subplots()
    ax.plot(range(len(kx_max)), kx_max, "-o")
    plt.show(block=False)

# Basic reconstruction
if traj_type == 0 and data_type == 'standard':
    image = mrr.recon_cart_rp(result)
elif data_type == 'standard':
    image = mrr.recon_noncart_rp(result)
else:
    image = 0
# Visualization
fig = mrr.view_image_montage(image, axis="x", title="MRI Reconstruction")


# -----------MRD CONVERT-----------#
if data_type == "gas exchange" or data_type == "calibration":
    pass
    # config_settings = {"ext_traj": True,
    #                    "trajorder": 1,
    #                    "contrast_order": [2, 1, 3],
    #                    "data_type": DataType.DIXON,
    #                    "xe_dissolved_offset_ppm": 202
    #                    }
    # XeGasExchange2XeCTCMRD.Gx2XeCTCMRD(
    #                     data_file=data_file,
    #                     raw_file=raw_file,
    #                     traj_file=traj_file,
    #                     config_settings=config_settings)
else:
    raw_data = p2m.Ph2Mrd(dlName=data_file, rlsName=raw_file)
    mrd_file, rlsData, dlData = raw_data.convert(".")

summary = mrc.summarize_ismrmrd_file(mrd_file)
written = mrc._write_json_report(summary, str(mrc._default_output_path(mrd_file)))


# -----------MRD RECONSTRUCT-----------#
# Read in data
mrd_dict = mrr.read_ismrmrd_file(mrd_file, lean=True)

# Basic data views
if mrd_dict["header"].encoding[0].trajectory.name == "CARTESIAN":
    readouts = np.squeeze(mrd_dict["data"])
    readouts = readouts.reshape(-1, readouts.shape[-1])
    fig, ax = plt.subplots()
    ax.plot(range(np.shape(readouts)[-1]), np.transpose(readouts, axes=(1,0)))
    plt.show(block=False)
else:
    k0 = mrd_dict["data"][..., 0]
    k0 = abs(k0.reshape(-1))
    fig, ax = plt.subplots()
    ax.plot(range(len(k0)), k0, "-o")
    plt.show(block=False)

    kx_max = mrd_dict["traj"][..., -1, 0]
    kx_max = abs(kx_max.reshape(-1))
    fig, ax = plt.subplots()
    ax.plot(range(len(kx_max)), kx_max, "-o")
    plt.show(block=False)

# Basic reconstruction
if mrd_dict["header"].encoding[0].trajectory.name == "CARTESIAN":
    image = mrr.recon_cart_mrd(mrd_dict)
else:
    image = mrr.recon_noncart_mrd(mrd_dict)

# Visualization
fig = mrr.view_image_montage(image, axis="x", title="MRI Reconstruction")


# file_a = "testdata\\IRC186-525_GX_MRD_Test\\xenon\\xenon_dixon_true.h5"
# file_b = "testdata\\IRC186-525_GX_MRD_Test\\xenon\\xenon_dixon.h5"

# report = mrc.compare_ismrmrd_files(path_a=file_a,
#                                    path_b=file_b,
#                                    compare_data=True)
# print(json.dumps(mrc._jsonable(report), indent=2, sort_keys=False))

# written = mrc._write_json_report(
#     report,
#     str(mrc._default_output_path(file_a, file_b)),
#     2,
# )
# print(f"JSON report written to: {written}", file=sys.stderr)
