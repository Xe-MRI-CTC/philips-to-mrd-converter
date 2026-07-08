import os
from pathlib import Path


def getTestData():
    # Set paths
    loc = Path(__file__).parent.absolute()
    data_path = os.path.join(loc, "testdata")

    # Set up test examples
    test_data = [
        {  # cartesian example
            "type": "standard",
            "name": "cartesian",
            "dl": os.path.join(data_path, "2DCartesian", "2DCartesian.data"),
            "rls": os.path.join(data_path, "2DCartesian", "2DCartesian.sin"),
            "traj": None
        },
        {  # spiral example
            "type": "standard",
            "name": "spiral",
            "dl": os.path.join(data_path, "2DSpiral", "2DSpiral.data"),
            "rls": os.path.join(data_path, "2DSpiral", "2DSpiral.sin"),
            "traj": None
        },
        {  # gas ex floret example
            "type": "gas exchange",
            "name": "gas exchange w/ floret",
            "dl": os.path.join(data_path, "3DFLORET_GXBonus", "3DFLORET_GXBonus.data"),
            "rls": os.path.join(data_path, "3DFLORET_GXBonus", "3DFLORET_GXBonus.sin"),
            "traj": os.path.join(data_path, "3DFLORET_GXBonus", "3DFLORET_GXBonus_Traj.sin")
        },
        {  # radial example
            "type": "standard",
            "name": "radial",
            "dl": os.path.join(data_path, "3DRadial", "3DRadial.data"),
            "rls": os.path.join(data_path, "3DRadial", "3DRadial.sin"),
            "traj": None
        },
        {  # gas ex w/ 2 echo example
            "type": "gas exchange",
            "name": "gas exchange w/ 2 echoes",
            "dl": os.path.join(data_path, "3DRadial_GX2Echo", "3DRadial_GX2Echo.data"),
            "rls": os.path.join(data_path, "3DRadial_GX2Echo", "3DRadial_GX2Echo.sin"),
            "traj": None
        },
        {  # gas ex w/ bonus spectra example
            "type": "gas exchange",
            "name": "gas exchange w/ bonus spectra",
            "dl": os.path.join(data_path, "3DRadial_GXBonus", "3DRadial_GXBonus.data"),
            "rls": os.path.join(data_path, "3DRadial_GXBonus", "3DRadial_GXBonus.sin"),
            "traj": os.path.join(data_path, "3DRadial_GXBonus", "3DRadial_GXBonus_Traj.sin")
        },
        {  # ctc gas ex example
            "type": "gas exchange",
            "name": "CTC gas exchange",
            "dl": os.path.join(data_path, "3DRadial_GXCTC", "3DRadial_GXCTC.data"),
            "rls": os.path.join(data_path, "3DRadial_GXCTC", "3DRadial_GXCTC.sin"),
            "traj": None
        },
        {  # XeCTC Calibration example
            "type": "calibration",
            "name": "CTC calibration",
            "dl": os.path.join(data_path, "XeCTC_Calibration", "XeCTC_Calibration.data"),
            "rls": os.path.join(data_path, "XeCTC_Calibration", "XeCTC_Calibration.sin"),
            "traj": None
        }
    ]

    return test_data
