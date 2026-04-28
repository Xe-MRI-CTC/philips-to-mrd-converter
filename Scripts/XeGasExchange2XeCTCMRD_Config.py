from enum import Enum, auto


class DataType(Enum):
    CALIBRATION = auto()
    DIXON = auto()
    UTE = auto()


class Config():
    '''Configuration class to correctly modify standard MRD conversion to XeCTC MRD'''

    def __init__(self):  
        '''
        Sets default values that are used as within the converter
        '''
        # default values hard coded for XeCTC acquisition at CCHMC

        # SET DEFAULT BASIC INFORMATION
        self.institution = 'CCHMC' # 'CCHMC' | 'Polarean'
        self.field_strength = 3.0
        self.H1resonanceFrequency_Hz = 127753955 
        self.orientation = 'Coronal'

        # SET DEFAULT TRAJECTORY INFORMATION
        self.ext_traj = False
        self.gr_delay = 1.25  # gradient delay used in calculating trajectories
        self.trajorder = 2  # trajectory ordering for radial acqusitions, 2 is haltoned spiral, 1 is 2D golden means, 0 is stock Philips

        # SET DEFAULT DATA TYPE INFORMATION
        self.data_type = DataType.CALIBRATION
        self.contrast_order = [1, 2]  # gas/diss
        self.flip_angle_gas = 0.5
        self.flip_angle_dis = 20.0
        self.cal_diss_acqs = 500  # all use 500 atm
        self.xe_dissolved_offset_ppm = 218.0

        # SET CUSTOM PARAMETERS
        self.prep_pulses = False
        self.gas_contam_removal = False  
        self.gas_contam_method = 'bonus_spec'
        self.exclude_bonus_spec = False  

        # SET DEBUG FLAG
        self.debug_mode = False

    def update(self, dl, rls, mrdHeader):
        '''
        Function to update the configuration parameters based on the data provided.
        Logic primarily based on scan name
        '''
        # Non scan-name-based updates:
        # UPDATE DATA TYPE
        # if multiple frequencies (stored as dynamics/repitition)
        if mrdHeader.encoding[0].encodingLimits.repetition.maximum > 0:
            self.data_type = DataType.DIXON     

        # Scan-name-based updates:
        if rls != None:
            scan_name = rls.header['sin']['scan_name'][0][0].lower()
        else:
            # data/list doesn't have scan name info so can't update based on name
            print('WARNING: No raw/lab/sin data prevents accurate update of parameters...')
            print('WARNING: Make sure to update config parameters manually.')
            return

        # UPDATE TRAJECTORY ORDER INFORMATION
        # two (V3) CPIR versions use golden means
        if 'Dissolved'.lower() in scan_name or \
           'CPIR'.lower() in scan_name:
            self.trajorder = 1
            self.gr_delay = +0.36
        if 'GX_CPIR_1pDIXON-V4'.lower() in scan_name:
            self.trajorder = 0
            self.gr_delay = +0.36 #likely not optimal

        # UPDATE DATA TYPE
        # If larger gamma indicating proton, set to UTE                              
        if float(rls.header['sin']['acq_gamma'][0][0]) > 42000.0:  # Proton
            self.data_type = DataType.UTE

        if self.data_type == DataType.UTE:
            return  # skip all xenon specific parameters

        # UPDATE INTERLEAVE ORDER
        # two CPIR versions collect diss/gas/off res
        if 'Dissolved'.lower() in scan_name or \
           'CPIR'.lower() in scan_name:
            self.contrast_order = [2, 1, 3]
        # FLORET version collect diss/gas
        if 'FLORET'.lower() in scan_name or \
           'GX_CPIR_1pDIXON-V4'.lower() in scan_name:
            self.contrast_order = [2, 1]

        # UPDATE FLIP ANGLE #can this be removed now is sin is accurate
        # Duke, FLORET, etc protocol uses smaller flip angle and corresponding TR in dixon
        if 'DukeIPF_Gas_Exchange'.lower() in scan_name or\
           'FLORET'.lower() in scan_name or \
           'Xenon_3D_radial_Dixon'.lower() in scan_name or \
           'GX_CPIR_1pDIXON-V4'.lower() in scan_name:
            self.flip_angle_dis = 15.0

        # Duke protocol sets dissolved between RBC and membrane for cal and dixon
        if 'Duke'.lower() in scan_name:
            self.xe_dissolved_offset_ppm = 208.0
        # two CPIR versions collect diss at 7143Hz
        if 'Dissolved'.lower() in scan_name or \
           'GX_CPIR_1pDIXON-V4'.lower() in scan_name:
            self.xe_dissolved_offset_ppm = 202.15
        # FLORET collect diss at 7143Hz
        if 'FLORET'.lower() in scan_name:
            self.xe_dissolved_offset_ppm = 202.15
            