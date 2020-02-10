#!/usr/bin/env python
# coding: utf-8

# In[ ]:

import numpy as np
import sys
from pathlib import Path
import glob
import os

# option to import from github folder
sys.path.insert(0, 'C:/Users/trose/Documents/GitHub/suite2p')
from suite2p import run_s2p

# set your options for running
ops = run_s2p.default_ops() # populates ops with the default options

#own imports
from ScanImageTiffReader import ScanImageTiffReader
from helpers import parse_SI_header as pSI #own
from shutil import copy, rmtree
from tqdm import tqdm

# ### File handling (PC / Mac)

# EXPERIMENT IDs
exp = ['62282', '62283', '62284', '62285', '62286', '62287', '62288', '62289', '62290', '62291', '62292', '62293', '62304', '62305', '62306', '62307', '62308', '62309', '62359', '62360', '62361', '62362', '62363', '62364']

if sys.platform == "win32":
    main_root = 'I:/David Laubender/Data/imaging data/DL_191024_6/ImagingData/' #location of original data
    ftemp = 'C:/temp/trose/suite2ptemp'        #fast disk
    # if len(exp)>1:
    #     froot = 'C:/temp/trose/s2p_tiff/'      #temp drive for concatenating exp tiffs
    #     fsave =  froot
    # else:
    froot = main_root                      #No concatenation of experiments -> can work on NW drive
    fsave = froot
elif sys.platform == "darwin":
    main_root = '/Volumes/archive_bonhoeffer_group$/David Laubender/Data/imaging data/DL_191024_6/ImagingData/' #location of original data
    ftemp = '/Users/trose/Data/temp/'          #fast disk
    # if len(exp)>1:
    #     froot = '/Users/trose/Data/s2p_tiff/'  #temp drive for concatenating exp tiffs
    #     fsave = froot
    # else:
    froot = main_root                      #No concatenation of experiments -> can work on NW drive
    fsave = froot
else:
    print('Unkknown Platform')

ftemp_tmp = ftemp + '/suite2p'
Path(ftemp_tmp).mkdir(parents=True, exist_ok=True)

for ids, val in enumerate(exp):

    files_fullpath = []
    files_name = []

    files = list(Path(froot).rglob('exp'+val+'*.tif')) #recursive

### continue here! extract path from files and updata froot and fsave_tmp

    fsave_tmp = froot + '/suite2p_exp' + exp[ids]
    Path(fsave_tmp).mkdir(parents=True, exist_ok=True)
    fsave = fsave_tmp

    print('- - - - - - - -')
    print('Processing exp' + val)
    print('- - - - - - - -')
    print('Folders used')
    print('Fast disk (ftemp) ' + ftemp)
    print('Data folder (froot) ' + froot)
    print('Save path (fsave) ' + fsave)
    print('- - - - - - - -')

    for n,f in enumerate(files):
        #files[n] = os.path.basename(f)
        files_fullpath.append(str(files[n]))
        files_name.append(f.name)

    files_fullpath = sorted(files_fullpath)
    files_name = sorted(files_name)



    with ScanImageTiffReader(files_fullpath[0]) as reader:
        header = (reader.description(0))
        mov_dim = (reader.shape())

    level = pSI.parse_SI_header_level(header)
    zoom = pSI.parse_SI_header_zoom(header)
    framerate = pSI.parse_SI_header_FrameRate(header)
    channels = pSI.parse_SI_header_Channels(header)
    volumes = pSI.parse_SI_header_Volumes(header)
    frames = pSI.parse_SI_header_Frames(header)
    frames_per_file = pSI.parse_SI_header_FramesPerFile(header)

    # account for multilevel acq where frames is 1
    if frames < volumes:
        frames = volumes


    bad_frames = np.linspace(0, 49, 50)
    np.save(fsave_tmp + '/bad_frames.npy', bad_frames)
    np.save(ftemp_tmp + '/bad_frames.npy', bad_frames)

    db = {
          'h5py': [], # a single h5 file path
          'h5py_key': 'data', # file paths
          'look_one_level_down': False, # whether to look in ALL subfolders when searching for tiffs
          'data_path': [froot],  # a list of folders with tiffs
                                 # (or folder of folders with tiffs if look_one_level_down is True, or subfolders is not empty)
          'delete_bin': False,
          'save_path0': fsave,
          'fast_disk': ftemp,
          'subfolders': [], # choose subfolders of 'data_path' to look in (optional)

          # main settings
          'nplanes': level, # each tiff has these many planes in sequence
          'nchannels': channels, # each tiff has these many channels per plane
          'functional_chan': 1, # this channel is used to extract functional ROIs (1-based)
          'tau':  1., # this is the main parameter for deconvolution
          'fs': framerate / level,  # sampling rate (PER PLANE - e.g. if you have 12 planes then this should be around 2.5)
          'preclassify': 0, # apply classifier before signal extraction with probability 0.5 (turn off with value 0)
          'frames_include:': -1,  #default: -1) if greater than zero, only frames_include frames are processed. useful for testing parameters on a subset of data.
          # output settings
          'save_mat': True, # whether to save output as matlab files
          'combined': True, # combine multiple planes into a single result /single canvas for GUI

          'num_workers': 0, # 0 to select num_cores, -1 to disable parallelism, N to enforce value
          'num_workers_roi': 0, # 0 to select number of planes, -1 to disable parallelism, N to enforce value
          'force_sktiff': False,
          # bidirectional phase offset
          'do_bidiphase': True,
          'bidiphase': 0,

          # registration settings
          'do_registration': 2, # whether to register data (2 forces re-registration)
          'two_step_registration:': True, #default: False) whether or not to run registration twice (for low SNR data). keep_movie_raw must be True for this to work.
          'keep_movie_raw': True,
          'nimg_init': 400, # subsampled frames for finding reference image
          'batch_size': 800, # number of frames per batch
          'maxregshift': 0.05, # max allowed registration shift, as a fraction of frame max(width and height)
          'align_by_chan' : 1, # when multi-channel, you can align by non-functional channel (1-based)
          'reg_tif': True, # whether to save registered tiffs
          'reg_tif_chan2': False, # whether to save channel 2 registered tiffs
          'subpixel' : 10, # precision of subpixel registration (1/subpixel steps)
          'smooth_sigma': 1.5, # ~1 good for 2P recordings, recommend >5 for 1P recordings
          'smooth_sigma_time': 2, # default: 0) standard deviation in time frames of the gaussian used to smooth the data before phase correlation is computed. Might need this to be set to 1 or 2 for low SNR data.
          'th_badframes': 2, # this parameter determines which frames to exclude when determining cropping - set it smaller to exclude more frames
          'pad_fft': False,

          # non rigid registration settings
          'nonrigid': True, # whether to use nonrigid registration
          'block_size': [128, 128], # block size to register (** keep this a multiple of 2 **)
          'snr_thresh': 1.5, # if any nonrigid block is below this threshold, it gets smoothed until above this threshold. 1.0 results in no smoothing
          'maxregshiftNR': 8, # maximum pixel shift allowed for nonrigid, relative to rigid

          # cell detection settings
          'roidetect': True, # whether or not to run ROI extraction
          'spatial_scale': 0, # 0: multi-scale; 1: 6 pixels, 2: 12 pixels, 3: 24 pixels, 4: 48 pixels
          'diameter': [9,12], #this is the main parameter for cell detection, 2-dimensional if Y and X are different (e.g. [6 12])
          'connected': False, # whether or not to keep ROIs fully connected (set to 0 for dendrites)
          'nbinned': 5000, # max number of binned frames for cell detection
          'max_iterations': 50, # maximum number of iterations to do cell detection
          'threshold_scaling': 2.0, # adjust the automatically determined threshold by this scalar multiplier
          'max_overlap': 0.9, # cells with more overlap than this get removed during triage, before refinement
          'high_pass': 100, # running mean subtraction with window of size 'high_pass' (use low values for 1P)
          'smooth_masks': True, # default: True) whether to smooth masks in final pass of cell detection. This is useful especially if you are in a high noise regime.

          # ROI extraction parameters
          'sparse_mode': False, #default: False) whether or not to use sparse_mode cell detection
          'inner_neuropil_radius': 2, # number of pixels to keep between ROI and neuropil donut
          'min_neuropil_pixels': 350, # minimum number of pixels in the neuropil
          'allow_overlap': False, # pixels that are overlapping are thrown out (False) or added to both ROIs (True)

          # channel 2 detection settings (stat[n]['chan2'], stat[n]['not_chan2'])
          'chan2_thres': 0.65, # minimum for detection of brightness on channel 2

          # deconvolution settings
          'baseline': 'maximin', # baselining mode (can also choose 'prctile')
          'win_baseline': 60., # window for maximin
          'sig_baseline': 10., # smoothing constant for gaussian filter
          'prctile_baseline': 8.,# optional (whether to use a percentile baseline)
          'neucoeff': .7,  # neuropil coefficient


          # List of tiffs to be loaded
          'tiff_list': files_name # list of tiffs in folder * data_path *!
    }

    opsEnd = run_s2p.run_s2p(ops=ops, db=db)

    rmtree(ftemp_tmp)
