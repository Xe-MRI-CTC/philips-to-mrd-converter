from enum import Enum, auto


class DataType(Enum):
    CALIBRATION = auto()
    DIXON = auto()
    UTE = auto()


class Config():
    '''Configuration class to correctly modify standard MRD conversion to XeCTC GasEx MRD'''

    def __init__(self):  # default values hard coded for XeCTC acquisition at CCHMC
        self.institution = 'CCHMC'
        self.field_strength = 3.0
        self.H1resonanceFrequency_Hz = 127753955
        self.orientation = 'Coronal'

        self.gr_delay = +2.5  # gradient delay used in calculating trajectories
        self.trajorder = 2  # trajectory ordering for radial acqusitions, 2 is halton randomized spiral, 1 is 2D golden means, 0 is stock Philips

        self.multi_echo = False
        self.ext_traj = False

        self.data_type = DataType.CALIBRATION
        self.contrast_order = [1, 2]  # gas/diss
        self.bonus_spec = False
        self.prep_pulses = False

        self.flip_angle_gas = 0.5
        self.flip_angle_dis = 20.0
        self.cal_diss_acqs = 500  # all use 500 atm
        self.xe_dissolved_offset_ppm = 218.0

    # default values hard coded for XeCTC acquisition at CCHMC
    def update(self, dl, rls, mrdHeader):
        # two CPIR versions use golden means
        if 'Dissolved'.lower() in rls.header['sin']['scan_name'][0][0].lower() or 'CPIR'.lower() in rls.header['sin']['scan_name'][0][0].lower():
            self.trajorder = 1
            self.gr_delay = +0.36

        # if multiple frequencies (stored as dynamics/repitition)
        if mrdHeader.encoding[0].encodingLimits.repetition.maximum > 0:
            self.data_type = DataType.DIXON
        if float(rls.header['sin']['acq_gamma'][0][0]) > 42000.0:  # Proton
            self.data_type = DataType.UTE

        # if multiple echoes stored as contrasts in mrd
        if mrdHeader.encoding[0].encodingLimits.contrast.maximum > 0:
            self.multi_echo = True

        if self.data_type == DataType.UTE:
            return  # skip all xenon specific parameters

        # 2 Mixes required to be cartesian; none of these should be so require an acquisition without 2 mixes (bonus spec)
        try:
            traj_type = int(rls.header['sin']['k_space_traj_type'][0][0])
        except:
            traj_type = 0
        if traj_type == 0 and not self.data_type == DataType.CALIBRATION:
            self.ext_traj = True

        # two CPIR versions collect diss/gas/off res
        if 'Dissolved'.lower() in rls.header['sin']['scan_name'][0][0].lower() or 'CPIR'.lower() in rls.header['sin']['scan_name'][0][0].lower():
            self.contrast_order = [2, 1, 3]
        # FLORET version collect diss/gas
        if 'FLORET'.lower() in rls.header['sin']['scan_name'][0][0].lower():
            self.contrast_order = [2, 1]

        # Mixes used to store bonus spec
        if int(rls.header['sin']['nr_mixes'][0][0]) > 1:
            self.bonus_spec = True

        # assume all trs the same and different frequencies are collected as different dynamics
        self.tr_factor = mrdHeader.encoding[0].encodingLimits.repetition.maximum + 1

        # Duke and FLORET protocol uses smaller flip angle and corresponding TR in dixon
        if 'DukeIPF_Gas_Exchange'.lower() in rls.header['sin']['scan_name'][0][0].lower():
            self.flip_angle_dis = 15.0
        if 'FLORET'.lower() in rls.header['sin']['scan_name'][0][0].lower():
            self.flip_angle_dis = 15.0

        # Duke protocol sets dissolved between RBC and membrane for cal and dixon
        if 'Duke'.lower() in rls.header['sin']['scan_name'][0][0].lower():
            self.xe_dissolved_offset_ppm = 208.0
        # two CPIR versions collect diss at 7143Hz
        if 'Dissolved'.lower() in rls.header['sin']['scan_name'][0][0].lower() or 'CPIR'.lower() in rls.header['sin']['scan_name'][0][0].lower():
            self.xe_dissolved_offset_ppm = 202.15
        # FLORET collect diss at 7143Hz
        if 'FLORET'.lower() in rls.header['sin']['scan_name'][0][0].lower():
            self.xe_dissolved_offset_ppm = 202.15
