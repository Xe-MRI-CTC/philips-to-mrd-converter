import numpy as np
from enum import Enum, auto

class DataType(Enum):
    CALIBRATION = auto()
    DIXON = auto()
    UTE = auto()

class Config():
    '''Configuration class to correctly modify standard MRD conversion to XeCTC GasEx MRD'''
    def __init__(self): # default values hard coded for XeCTC acquisition at CCHMC
        self.institution = 'CCHMC'
        self.field_strength = 3.0
        self.H1resonanceFrequency_Hz = 127753955
        self.orientation = 'Coronal'

        self.gr_delay = +2.5 # gradient delay used in calculating trajectories
        self.traj_order = 2 # trajectory ordering for radial acqusitions, 2 is halton randomized spiral, 1 is 2D golden means, 0 is stock Philips

        self.multi_echo = False
        self.ext_traj = False
        
        self.data_type = DataType.CALIBRATION
        self.contrast_order = [1, 2] # gas/diss
        self.bonus_spec = False

        self.flip_angle_gas = 0.5
        self.flip_angle_dis = 20.0
        self.cal_diss_acqs = 500 # all use 500 atm
        self.xe_dissolved_offset_ppm = 218.0

    def update(self, dl, rls, mrdHeader): # default values hard coded for XeCTC acquisition at CCHMC
        if 'Dissolved'.lower() in rls.header['sin']['scan_name'][0][0].lower() or 'CPIR'.lower() in rls.header['sin']['scan_name'][0][0].lower(): # two CPIR versions use golden means
            self.traj_order = 1
            self.gr_delay = +0.36

        if  mrdHeader.encoding[0].encodingLimits.repetition.maximum > 0: # if multiple frequencies (stored as dynamics/repitition)
            self.data_type = DataType.DIXON
        if float(rls.header['sin']['acq_gamma'][0][0]) > 42000.0: # Proton
            self.data_type = DataType.UTE

        if  mrdHeader.encoding[0].encodingLimits.contrast.maximum > 0: # if multiple echoes stored as contrasts in mrd
            self.multi_echo = True
        
        if self.data_type == DataType.UTE: 
            return # skip all xenon specific parameters

        if 'Dissolved'.lower() in rls.header['sin']['scan_name'][0][0].lower() or 'CPIR'.lower() in rls.header['sin']['scan_name'][0][0].lower(): # two CPIR versions are cartesian (no trajs) for added spec
            self.ext_traj = True

        if 'Dissolved'.lower() in rls.header['sin']['scan_name'][0][0].lower() or 'CPIR'.lower() in rls.header['sin']['scan_name'][0][0].lower(): # two CPIR versions collect diss/gas/off res
            self.contrast_order = [2, 1, 3]

        if 'Dissolved'.lower() in rls.header['sin']['scan_name'][0][0].lower() or 'CPIR'.lower() in rls.header['sin']['scan_name'][0][0].lower(): # two CPIR versions containing bonus spectra
            self.bonus_spec = True

        self.tr_factor = mrdHeader.encoding[0].encodingLimits.repetition.maximum + 1 # assume all trs the same and different frequencies are collected as different dynamics

        if 'DukeIPF_Gas_Exchange'.lower() in rls.header['sin']['scan_name'][0][0].lower(): # Duke protocol uses smaller flip angle and corresponding TR in dixon
            self.flip_angle_dis = 15.0 

        if 'Duke'.lower() in rls.header['sin']['scan_name'][0][0].lower(): # Duke protocol sets dissolved between RBC and membrane for cal and dixon
            self.xe_dissolved_offset_ppm = 208.0
        if 'Dissolved'.lower() in rls.header['sin']['scan_name'][0][0].lower() or 'CPIR'.lower() in rls.header['sin']['scan_name'][0][0].lower(): # two CPIR versions collect diss at 7143Hz
            self.xe_dissolved_offset_ppm = 202.15






 
