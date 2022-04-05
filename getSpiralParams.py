#!/usr/bin/env python
"""This is a library of routines for reading in and
populating a dictionary of the parameters needed for
spiral reconstruction.
"""

import os
import sys
import numpy as np
import re
from collections import OrderedDict

class Scanner(object):
    def __init__(self,ver):
        self.ver = ver

# match the correct filename case
#   * the existence check is case insensitive
#   * need the exact filename for opening a file
#   * the basename is chosen by the user so we just need to check for ext case
#   * this should be used when the file extension is guessed
def filename_extcase(fn):
    pn, ext = os.path.splitext(fn)
    bn = os.path.basename(pn)
    if bn+ext.lower() in os.listdir(os.path.dirname(fn)):
        return pn+ext.lower()
    elif bn+ext.upper() in os.listdir(os.path.dirname(fn)):
        return pn+ext.upper()
    return ''


# BNI spiral specific code for parsing parameter .txt files
# Files written in the YAC format are parsed and put into a dict.
#   See mrmethods/methpdf/src/mparc.cpp for details.
#
# A class to aid in writing YAC Ain't CSV (YAC) files.
#
#  File Format:
#
#  <param-name>: <value1>, <value2>, ... <value_n>
#                                                 ^ NO trailing commas!
#
#  NOTE: 1) Comma separated value lists can also convey the number of values
#           without having to pass an added 'length' parameter.
#
#        2) The first colon separates the param-name from the values,
#           additional colons will be treated as text.
#
#        3) All the values will eventually be converted to a float header.
#           Since text values will eventually be enumerated, numeric flags are
#           easier to pass at this point.
#
#        4) All comment lines need to start with a colon ":".
#
#        5) No blank lines.
#
# DICT OUTPUT:
#
#   The output dictionary is effectively a direct translation of the YAC file.
#   The <param-name>s are converted to keys and the <valuen>s (even single
#   values) are put into list objects.  The values are all kept as strings, so
#   the user must recast them as necessary. The dict is one-level and the value
#   lists are one-level.
def readParms(filename):
    '''
        Parse a *.txt (or *.yac) header accompanying the BNI spiral *.data
        *.list files. Written in the YAC file format.
    '''

    # open the txt file
    if not type(filename) is str:
        print("Input filename is not a string.")
        sys.exit(1)

    if str(os.path.splitext(filename)[1]).lower() not in ['.txt', '.yac']:
        print("input filename is not a .txt or .yac file")
        sys.exit(1)

    # opens the file
    try:
        fil = open(filename_extcase(filename), 'r', encoding='ascii')
    except IOError:
        print('cannot open .txt (or .yac) file ', filename)
        sys.exit(1)

    # read in all the lines at one time for processing
    lines = fil.readlines()
    fil.close()

    # remove white space and delimiters etc..
    # whitespace
    lines = [re.sub('\s', '', line) for line in lines]
    lines = [line for line in lines if line.strip() != '']
    lines = [re.split(':', line, maxsplit=1) for line in lines]

    # create the header dictionary
    hdr = dict()
    for line in lines:
        hdr[line[0]] = line[1]

    # split dictionary values into sub arrays
    # NOTE: this will put single values into a dict as well
    for key, value in hdr.items():
        hdr[key] = value.split(',')

    # so that GPI nodes can verify where this dict came from
    hdr['headerType'] = 'BNIspiral'

    # remove comment entries
    del hdr['']

    return hdr


def readAdhoc(adhoc_float,adhoc_int,scanner):
    '''
        Parse the SpiralAdhocR56.h header and create a dictionary from the .sin
        file ad_hoc_float_array values.
    '''

    # opens the file if it exists
    filename = os.path.dirname(__file__)+'/SpiralAdhoc'+scanner.ver+'.h'
    try:
        fil = open(filename, 'r')
    except IOError:
        print('cannot open', filename)
        sys.exit(1)

    # read in all the lines at one time for processing
    lines = fil.readlines()
    fil.close()

    # remove lines that aren't part of AD HOC FLOAT or AD HOC INT
    begin = [num for num, line in enumerate(lines, 1)
             if(line.find('AD HOC FLOAT') != -1)][0] - 1
    begin_int = [num for num, line in enumerate(lines, 1)
             if(line.find('AD HOC INT') != -1)][0] - 1         
    end = [num for num, line in enumerate(lines, 1)
           if(line.find('#endif') != -1)][0] - 1
    lines_float = lines[begin:begin_int]
    lines_int = lines[begin_int:end]

    # remove lines w/o #define
    lines_float = [line for line in lines_float if(line.find('#define') != -1)]
    lines_int = [line for line in lines_int if(line.find('#define') != -1)]

    # split on spaces, remove #define
    lines_float = [re.split('\s+', line.strip(
        '#define\r\n').strip()) for line in lines_float]
    lines_int = [re.split('\s+', line.strip(
        '#define\r\n').strip()) for line in lines_int]

    # create dictionary with labels from SpiralAdhocR56.h and values from
    # the ad_hoc_float_array
    params = OrderedDict()
    for line in lines_float:
        params[line[0]] = adhoc_float[int(line[1])]
    for line in lines_int:
        params[line[0]] = adhoc_int[int(line[1])]


    return params


def readSpiralCoords(filename, header, scanner):
    filename_raw = filename+'.raw'

    # opens the file
    try:
        fil_raw = open(filename_extcase(filename_raw), 'rb')
    except IOError:
        print('cannot open', filename_raw)

    # skip 512 bytes
    fil_raw.seek(512)

    # label file block size
    lab_pos = 0  # seek position

    # find indices and seek offsets for coord data
    lab = header['lab']
    std_case = lab['label_type'] == b'LABEL_TYPE_STANDARD'
    if scanner.ver == 'R56':
        coord_case = lab['control'] == b'CTRL_TRAJ_DATA'
    else:
        coord_case = lab['control'] == b'CTRL_RESUME'
    coord_ind = (std_case & coord_case).nonzero()[0]
        
    # print("coord_ind:", np.sum(coord_ind))
    print("scanner version:", scanner.ver)

    # find seek positions in raw file
    sin = header['sin']
    isMira = sin['isMira']
    format_vals = lab['raw_format'] == 4
    format_vals |= lab['raw_format'] == 6
    if isMira:
        seek_array = np.cumsum(np.where(format_vals, lab['coded_data_size'],
                               lab['data_size']))
    else:
        seek_array = np.cumsum(lab['data_size'])
    seek_array = np.insert(seek_array, 0, 0) + 512
    seek_vals = seek_array[coord_ind]

    # set up coordinate array
    nsamp = (int(sin['non_cart_max_encoding_nrs'][0][0]) -
             int(sin['non_cart_min_encoding_nrs'][0][0]) + 1)
    ncoord = len(coord_ind)
    if(scanner.ver == 'R56'):
        ncoord *= 3
    if ncoord > 0:
        coord_array = np.zeros([ncoord, nsamp])
    else:
        coord_array = None

    cur_coord = 0
    # read the coordinates from the raw file
    for seek_ind in range(len(seek_vals)):
        seek_offset = int(seek_vals[seek_ind])
        lab_pos = coord_ind[seek_ind]
        act_data_size = int(lab['data_size'][lab_pos])

        fil_raw.seek(seek_offset)
        bytebuff = fil_raw.read(act_data_size)
        temp_data = np.fromstring(bytebuff, dtype=np.float32)
        if (scanner.ver == 'R56'):
            cpxData = temp_data[0:3*nsamp]
        else:
            cpxData = temp_data[0:nsamp]
            coord_array[cur_coord,:] = cpxData
            cur_coord = cur_coord + 1

    if(scanner.ver == 'R56'):
        coord_array = cpxData.reshape((nsamp,ncoord)).T
    # close the raw binary file
    fil_raw.close()
    return(coord_array)


def processSpiralParams(header, filename, scanner):
    spparams = dict()
    spparams['headerType'] = 'spparams'

    hdr = None
    if header['headerType'] == 'BNIspiral':
        hdr = header
    elif 'BNIspiral' in header:
        hdr = header['BNIspiral']

    sin = None
    spiral_traj = False
    if header['headerType'] == 'lab-sin':
        sin = header['sin']
        if 'k_space_traj_type' in sin:
            spiral_traj = (int(float(sin['k_space_traj_type'][0][0])) == 2)

    if hdr is not None:
        # adhoc parameters sent to recon
        spparams['NR_SPIRAL_ARMS'] = float(hdr['spARMS'][0])
        spparams['SPIRAL_TYPE'] = int(float(hdr['spSTYPE'][0]))
        spparams['SPIRAL_IN_OUT_TYPE'] = int(float(hdr['spINOUT_ON'][0]))
        spparams['SPIRAL_INOUT_OPT'] = int(float(hdr['spINOUT_OPT'][0]))
        spparams['DWELL_US'] = 1000. * float(hdr['spDWELL'][0])
        spparams['NR_GRD_SAMPLES_RAMP_DOWN'] =\
            (int(float(hdr['spgrad_nb'][0]))
             - int(float(hdr['spgrad_na'][0])))
        if 'spEXTRA_GRAD_PNT' in hdr:
            spparams['ZERO_GRD_POINTS'] = int(float(
                hdr['spEXTRA_GRAD_PNT'][0]))

        # f0 navigator values
        if 'spF0nav_LONG' in hdr:
            f0ref = float(hdr['spF0nav_LONG'][0])
            df01 = float(hdr['spF0nav_LONG'][1]) - f0ref
            df02 = float(hdr['spF0nav_LONG'][2]) - f0ref
            # if an f0 value is 0, df0 will be > 1000
            if abs(df01) < 1000. and abs(df02) < 1000.:
                spparams['F0_OFFSET'] = df01
                spparams['F0_OFFSET2'] = df02
            elif abs(df01) < 1000.:
                spparams['F0_OFFSET'] = df01
            elif abs(df02) < 1000.:
                spparams['F0_OFFSET'] = df02
            else:
                spparams['F0_OFFSET'] = 0.
        else:
            spparams['F0_OFFSET'] = 0.

        # parameters for generating the base spiral gradient waveform
        spparams['MAX_SLEW'] = float(hdr['spSLEWMAX'][0])
        spparams['MAX_GRAD'] = float(hdr['spGMAX'][0])
        spparams['GAMMA'] = float(hdr['spGAMMA'][0])
        fovxy = 100. * float(hdr['spFOVXY'][0])
        print("fovxy =", fovxy)
        fovz = 100. * float(hdr['spFOVZ'][0])
        spparams['FOV_CM'] = np.array([fovxy, fovxy, fovz])
        resxy = 100. * float(hdr['spRESXY'][0])
        print("resxy =", resxy)
        resz = 100. * float(hdr['spRESZ'][0])
        spparams['RES_CM'] = np.array([resxy, resxy, resz])
        if 'trures_fac_xy' in hdr:
            spparams['TRUE_RESXY'] = float(hdr['trures_fac_xy'][0])
        if 'trures_fac_z' in hdr:
            spparams['TRUE_RESZ'] = float(hdr['trures_fac_z'][0])
        spparams['UNDERSAMP_TYPE'] = int(float(hdr['spUSTYPE'][0]))
        if 'spUSR0' in hdr:
            spparams['MIN_UNDERSAMP'] = float(hdr['spUSR0'][0])
        spparams['MAX_UNDERSAMP'] = float(hdr['spUSR'][0])
        spparams['UNDERSAMP_START'] = float(hdr['spUS0'][0])
        spparams['UNDERSAMP_END'] = float(hdr['spUS1'][0])
        spparams['NON_SPIRAL_GRAD_TYPE'] = int(float(hdr['spGTYPE'][0]))
        if 'spMGFRQ' in hdr:
            spparams['MAX_GRAD_FREQ_kHZ'] = float(hdr['spMGFRQ'][0])
        if 'spSLWIN' in hdr:
            spparams['GRAD_SLEW_WIN_US'] = float(hdr['spSLWIN'][0])
        spparams['T2_MATCH_MS'] = float(hdr['spT2MATCH'][0])

        # off-center FOV offsets
        m_offc = (.125 * float(hdr['m_offc'][0]) / spparams['RES_CM'][0])
        p_offc = (.125 * float(hdr['p_offc'][0]) / spparams['RES_CM'][1])
        if spparams['RES_CM'][2] > 0.:
            s_offc = (.1 * float(hdr['s_offc'][0]) / spparams['RES_CM'][2])
        else:
            s_offc = 0.
        spparams['FOV_OFFC_PIXELS'] = np.array([m_offc, p_offc, s_offc])

        # additional spiral parameters
        spparams['SPIRAL_SAMP'] = int(float(hdr['spREADPTS'][0]))
        if 'spf_NSA' in hdr:
            spparams['FRACTIONAL_NSA'] = float(hdr['spf_NSA'][0])

        # FLORET parameters
        if 'spHUBS' in hdr:
            spparams['FLORET_HUBS'] = int(float(hdr['spHUBS'][0]))
        if 'spALPHA0' in hdr:
            spparams['FLORET_ALPHA0'] = float(hdr['spALPHA0'][0])
        if 'spFLORETbin' in hdr:
            spparams['FLORET_REBIN'] = int(float(hdr['spFLORETbin'][0]))

        # calculate the echo-times
        if 'gnNUMECHOES' in hdr:
            necho = int(float(hdr['gnNUMECHOES'][0]))
            spparams['NUM_ECHOES'] = necho
        else:
            necho = 1
        if 'spTESHIFT_ON' in hdr:
            teshift_on = int(float(hdr['spTESHIFT_ON'][0]))
        if 'spTESHIFT' in hdr and teshift_on == 1:
            te_shifts = np.array(hdr['spTESHIFT'], dtype=np.float32)
        else:
            te_shifts = np.zeros(1)
        if 'gnTE' in hdr:
            base_te = np.array(hdr['gnTE'][0:necho], dtype=np.float32)
        else:
            base_te = np.zeros(1)
        echo_time = base_te + te_shifts
        spparams['TE'] = echo_time

        # other scan info
        if 'gnTR' in hdr:
            spparams['TR'] = float(hdr['gnTR'][0])
        if 'gnNUMSLICE' in hdr:
            spparams['NUM_SLICES'] = int(float(hdr['gnNUMSLICE'][0]))
        if 'spTSESPIRAL_ETL' in hdr:
            spparams['spTSESPIRAL_ETL'] = int(float(hdr['spTSESPIRAL_ETL'][0]))

    if header['headerType'] == 'lab-sin' and spiral_traj:
        # get the exported coordinates
        # scanner = Scanner('R56')
        FLORETR56 = False # need to default this in case of data with no ad_hoc values
        spparams['COORDS'] = readSpiralCoords(filename, header, scanner)

        # adhoc parameters sent to recon (read appropriate SpiralAdhoc .h file if present)
        defpath = os.path.dirname(__file__)+'/SpiralAdhoc'+scanner.ver+'.h'
        if os.path.isfile(defpath) and ('ad_hoc_float_array' in sin or 'ad_hoc_int_array' in sin):
            # find indices where ad_hoc array has non-zero values
            # copy existing adhoc float values to an array
            adhoc_float = np.zeros(24)
            if 'ad_hoc_float_array' in sin:
                float_indices = np.array((sin['ad_hoc_float_array'][1]),
                                   dtype=int)[:, 0] - 1
                adhoc_float[float_indices] = np.array(sin['ad_hoc_float_array'][0],
                                              dtype=float)
                                      
            # repeat for ad_hoc_int_array, if present
            adhoc_int = np.zeros(24,dtype=int)
            if 'ad_hoc_int_array' in sin:
                int_indices = np.array((sin['ad_hoc_int_array'][1]),
                                   dtype=int)[:, 0] - 1
                adhoc_int[int_indices] = np.array(sin['ad_hoc_int_array'][0],
                                          dtype=int)
                                      
            # create parameter dictionary
            spparams.update(readAdhoc(adhoc_float,adhoc_int,scanner))
            try: #R56 FLORET
                if int(sin['spiral_trajectory_shape'][0][0]) == 3 and (scanner.ver == 'R56'):
                    #overwrite adhoc params not found in FLORET patch
                    spparams['SPIRAL_TYPE'] = 3 #FLORET
                    spparams.pop('NR_SPIRAL_ARMS', None)
                    FLORETR56 = True
                else:
                    FLORETR56 = False
            except:
                FLORETR56 = False
        else:
            print("Warning! No ad_hoc values in sin file. The necessary info"
                  " may not be available for spiral recon.")

        # resolution, fov, and off-center FOV
        recon_mtx = np.array(sin['recon_resolutions'][0][0:3], dtype=np.int32)
        voxel_size = np.array(sin['voxel_sizes'][0], dtype=float)
        if FLORETR56: #R56 FLORET - Force Isotropic
            recon_mtx[2] = recon_mtx[0]

        if(scanner.ver == 'R56'):
            fov = (voxel_size * recon_mtx).astype(np.float32)
            if FLORETR56: #R56 FLORET - Force Isotropic
                fov[2] = fov[0]
        else:
            fov = (voxel_size * recon_mtx).astype(np.int32)
        spparams['FOV_CM'] = fov / 10.

        oversamp = np.array(sin['oversample_factors'][0], dtype=float)
        if(scanner.ver == 'R56'):
            spparams['OVER_SAMP'] = oversamp 
            encoding_numbers =\
                (np.array(sin['max_encoding_numbers'][0], dtype=np.int32) -
                  np.array(sin['min_encoding_numbers'][0], dtype=np.int32) + 1)[0:3]
            res = np.where(encoding_numbers != 0., fov[0:3]*oversamp[0:3]/ encoding_numbers, 0.)
            if FLORETR56:
                res[2] = res[0]
        else:
            encoding_numbers =\
                ((np.array(sin['max_encoding_numbers'][0], dtype=np.int32) -
                  np.array(sin['min_encoding_numbers'][0], dtype=np.int32) + 1) /
                  oversamp)[0:3]
            res = np.where(encoding_numbers != 0., fov / encoding_numbers, 0.)

        offc = np.array(sin['location_center_coordinates'][0][0:3],
                        dtype=float)
        omatrix = np.array(sin['location_matrices'][0][0:9], dtype=float)
        omatrix.shape = [3, 3]
        omatrix = np.ascontiguousarray(np.transpose(omatrix))
        omatrix.shape = [9]

        m_offc = (offc[0]*omatrix[0] + offc[1]*omatrix[1]
                  + offc[2]*omatrix[2]) / res[0]
        p_offc = (offc[0]*omatrix[3] + offc[1]*omatrix[4]
                  + offc[2]*omatrix[5]) / res[1]
        s_offc = np.where(res[2] > 0., (offc[0]*omatrix[6] + offc[1]*omatrix[7]
                  + offc[2]*omatrix[8]) / res[2], 0.)
        if FLORETR56:#not robust by any means but a start until more data
            #not a huge deal as offc van be changed via widget
            m_offc = (offc[1]*-1)/(voxel_size[0])
            p_offc = (offc[2]*1)/(voxel_size[1])
            s_offc = (offc[0]*-1)/(voxel_size[2]*1.25)

        spparams['FOV_OFFC_PIXELS'] = np.array([m_offc, p_offc, s_offc])
        spparams['O_MATRIX'] = omatrix
        if(scanner.ver == 'R56'):
            if FLORETR56: #R56 FLORET
                res[0:3] = res[0:3]*1.25
            else:
                res[0:2] = res[0:2]*1.25
        else:
            res[0:2] *= 1.25
            if 'SPIRAL_TYPE' in spparams and spparams['SPIRAL_TYPE'] > 1:
                res[2] *= 1.25
        spparams['RES_CM'] = res/ 10 

        # calculate the echo-times
        necho = int(float(sin['nr_echoes'][0][0]))
        spparams['NUM_ECHOES'] = necho
        num_shifts = int(float(sin['nr_rows'][0][0]))
        if int(sin['sequence_types'][0][0]) == 0: #Spin Echo
            base_te = np.zeros(necho)
        else: #FFE
            base_te = np.array(sin['echo_times'][0][0:necho], dtype=np.float32)

        if 'dixon_enable' in sin: #dixon
            dixon_offset_te = float(sin['dixon_first_echo_time'][0][0])
            teshifts = (float(sin['dixon_delta_te'][0][0]) *
                        np.arange(num_shifts))
            echo_time = base_te + teshifts
        else: #non-dixon
            echo_time = base_te
        spparams['TE'] = echo_time

        # additional spiral parameters
        if(scanner.ver == 'R56'):
            spparams['SPIRAL_IN_OUT_TYPE'] = int(sin['spiral_in_out_type'][0][0]) 
            #spparams['X_DELAY_US'] = float(sin['gradient_delays'][0][0])
            #spparams['Y_DELAY_US'] = float(sin['gradient_delays'][0][1])
            #spparams['Z_DELAY_US'] = float(sin['gradient_delays'][0][2])
            spparams['DWELL_US'] = float(sin['sample_time_interval'][0][0])
            spparams['NR_GRD_SAMPLES_RAMP_DOWN'] = int(sin['spiral_nr_grd_smpls_rmp_dn'][0][0])
            spparams['ZERO_GRD_POINTS'] = int(sin['spiral_inout_zero_grd_pts'][0][0])
        nsamp = (int(sin['non_cart_max_encoding_nrs'][0][0]) -
                 int(sin['non_cart_min_encoding_nrs'][0][0]) + 1)
        spparams['SPIRAL_SAMP'] = nsamp
        if 'NR_SPIRAL_ARMS' in spparams:
            fNSA = ((int(sin['non_cart_max_encoding_nrs'][0][1]) -
                     int(sin['non_cart_min_encoding_nrs'][0][1]) + 1) /
                    spparams['NR_SPIRAL_ARMS'])
            spparams['FRACTIONAL_NSA'] = fNSA
        else:
            spparams['NR_SPIRAL_ARMS'] = (int(sin['non_cart_max_encoding_nrs'][0][1]) -
                     int(sin['non_cart_min_encoding_nrs'][0][1]) + 1) 

        # other scan info
        spparams['GAMMA'] = float(sin['acq_gamma'][0][0]) / 1000.
        if 'repetition_times' in sin:
            spparams['TR'] = float(sin['repetition_times'][0][0])
        if 'nr_locations' in sin:
            spparams['NUM_SLICES'] = int(float(sin['nr_locations'][0][0]))
        if FLORETR56: #R56 FLORET
            spparams['NUM_SLICES'] = int(sin['scan_resolutions'][0][0])
    
    return(spparams)


# for testing functions individually
if __name__ == '__main__':
    spparms = readParms('./spiraltxt.yac')
    for k, v in spparms.items():
        print("\'"+k+"\'", v)
