#!/usr/bin/env python
# coding: utf-8

# Tobias Rose 2020

import numpy as np
from pathlib import Path
import glob
import os
import sys
from suite2p import run_s2p
ops = run_s2p.default_ops() # populates ops with the default options

# option to import from github folder
sys.path.insert(0, 'C:/Users/trose/Documents/GitHub/suite2p')
from shutil import copy, rmtree
from tqdm import tqdm
from helpers import parse_SI_header as pSI
from ScanImageTiffReader import ScanImageTiffReader


def batch_s2p( runops ):
    "Batch run for s2p'ing through experiments. Settings are currently for axons"

    # populates s2p ops with the default options

    # Unpack option list
    exp = runops.get('exp')   #shorter call is runops['exp'], of course
    concat = runops.get('concat')
    main_root = runops.get('main_root')

    adata = runops.get('adata')
    ftiff = runops.get('ftiff')
    ftemp = runops.get('ftemp')

    darkframes = runops.get('darkframes')
    readfiles = runops.get('readfiles')

    # make subfolders
    ftemp_tmp = os.path.join(ftemp, 'suite2p')
    Path(ftemp_tmp).mkdir(parents=True, exist_ok=True)
    print('- - - - - - - -')
    print('Made directory (ftemp_tmp) ' + ftemp_tmp)

    # Concatenation case folder settings -> s2p root folder is the concat tiff folder
    if concat:
        froot = ftiff
        fsave_tmp =  os.path.join(froot, 'suite2p_exp')
        fsave = fsave_tmp
        Path(fsave_tmp).mkdir(parents=True, exist_ok=True)
        print('Made directory (fsave_tmp) ' + fsave_tmp)
        print('- - - - - - - -')
        print('Concatenating all experiments in single tiff folder')

        bad_frames = np.linspace(0, darkframes-1, darkframes)
        np.save(fsave_tmp + '/bad_frames.npy', bad_frames)
        np.save(ftemp_tmp + '/bad_frames.npy', bad_frames)

        files_fullpath = []
        files_name = []
        for val in exp:
            files = list(Path(main_root).rglob('exp'+val+'*.tif')) #recursive search over main_root
            files = sorted(files)
            files = files[0:readfiles]
            for n,f in enumerate(files):
                copy(files[n], froot)
                targetstring = [str(files[n]), str(froot)]
                print('copied {} to {}'.format(*targetstring))
                files_fullpath.append(str(files[n]))
                files_name.append(f.name)
        files_fullpath = sorted(files_fullpath)
        files_name = sorted(files_name)

        info = si_info( files_fullpath )

        runops = {
            'exp': exp,
            'files_name': files_name,
            'main_root': main_root,
            'adata': adata,
            'ftemp': ftemp,
            'ftiff': ftiff,
            'froot': froot,
            'fsave': fsave,
            'level': info.get('level'),
            'channels': info.get('channels'),
            'framerate': info.get('framerate'),
            }
        db = s2p_ownops(runops)
        opsEnd = run_s2p.run_s2p(ops=ops, db=db)
    else:
        for ids, val in tqdm(enumerate(exp)):

            files_fullpath = []
            files_name = []

            files = list(Path(main_root).rglob('exp'+val+'*.tif')) #recursive
            files = files[0:readfiles]

            froot = os.path.join(os.path.dirname(files[0]))

            fsave_tmp = os.path.join(froot, 'suite2p_exp' + exp[ids])
            fsave = fsave_tmp

            Path(fsave_tmp).mkdir(parents=True, exist_ok=True)
            print('Made directory (fsave_tmp) ' + fsave_tmp)

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

            bad_frames = np.linspace(0, darkframes-1, darkframes)
            np.save(fsave_tmp + '/bad_frames.npy', bad_frames)
            np.save(ftemp_tmp + '/bad_frames.npy', bad_frames)

            info = si_info( files_fullpath )

            runops = {
                'exp': exp,
                'files_name': files_name,
                'main_root': main_root,
                'adata': adata,
                'ftemp': ftemp,
                'ftiff': ftiff,
                'froot': froot,
                'fsave': fsave,
                'level': info.get('level'),
                'channels': info.get('channels'),
                'framerate': info.get('framerate'),
                }

            db = s2p_ownops(runops)
            opsEnd = run_s2p.run_s2p(ops=ops, db=db)

def batch_s2p_ROI( runops ):
    "Batch run for s2p'ing through experiments. Settings are currently for axons"

    # populates s2p ops with the default options

    # Unpack option list
    exp = runops.get('exp')   #shorter call is runops['exp'], of course
    concat = runops.get('concat')
    main_root = runops.get('main_root')

    adata = runops.get('adata')
    ftiff = runops.get('ftiff')
    ftemp = runops.get('ftemp')

    darkframes = runops.get('darkframes')
    readfiles = runops.get('readfiles')

    # make subfolders
    ftemp_tmp = os.path.join(ftemp, 'suite2p')
    Path(ftemp_tmp).mkdir(parents=True, exist_ok=True)
    print('- - - - - - - -')
    print('Made directory (ftemp_tmp) ' + ftemp_tmp)

    # Concatenation case folder settings -> s2p root folder is the concat tiff folder
    if concat:
        froot = ftiff
        fsave_tmp =  os.path.join(froot, 'suite2p_exp')
        fsave = fsave_tmp
        Path(fsave_tmp).mkdir(parents=True, exist_ok=True)
        print('Made directory (fsave_tmp) ' + fsave_tmp)
        print('- - - - - - - -')
        print('Concatenating all experiments in single tiff folder')

        bad_frames = np.linspace(0, darkframes-1, darkframes)
        np.save(fsave_tmp + '/bad_frames.npy', bad_frames)
        np.save(ftemp_tmp + '/bad_frames.npy', bad_frames)

        files_fullpath = []
        files_name = []
        for val in exp:
            files = list(Path(main_root).rglob('exp'+val+'*.tif')) #recursive search over main_root
            files = sorted(files)
            files = files[0:readfiles]
            for n,f in enumerate(files):
                copy(files[n], froot)
                targetstring = [str(files[n]), str(froot)]
                print('copied {} to {}'.format(*targetstring))
                files_fullpath.append(str(files[n]))
                files_name.append(f.name)
        files_fullpath = sorted(files_fullpath)
        files_name = sorted(files_name)

        info = si_info( files_fullpath )

        runops = {
            'exp': exp,
            'files_name': files_name,
            'main_root': main_root,
            'adata': adata,
            'ftemp': ftemp,
            'ftiff': ftiff,
            'froot': froot,
            'fsave': fsave,
            'level': info.get('level'),
            'channels': info.get('channels'),
            'framerate': info.get('framerate'),
            }
        db = s2p_ownops(runops)
        opsEnd = run_s2p.run_s2p(ops=ops, db=db)
    else:
        for ids, val in tqdm(enumerate(exp)):

            files_fullpath = []
            files_name = []

            files = list(Path(main_root).rglob('exp'+val+'*.tif')) #recursive
            files = files[0:readfiles]

            froot = os.path.join(os.path.dirname(files[0]))

            fsave_tmp = os.path.join(froot, 'suite2p_exp' + exp[ids])
            fsave = fsave_tmp

            Path(fsave_tmp).mkdir(parents=True, exist_ok=True)
            print('Made directory (fsave_tmp) ' + fsave_tmp)

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

            bad_frames = np.linspace(0, darkframes-1, darkframes)
            np.save(fsave_tmp + '/bad_frames.npy', bad_frames)
            np.save(ftemp_tmp + '/bad_frames.npy', bad_frames)

            info = si_info( files_fullpath )

            runops = {
                'exp': exp,
                'files_name': files_name,
                'main_root': main_root,
                'adata': adata,
                'ftemp': ftemp,
                'ftiff': ftiff,
                'froot': froot,
                'fsave': fsave,
                'level': info.get('level'),
                'channels': info.get('channels'),
                'framerate': info.get('framerate'),
                }

            db = s2p_ownops(runops)
            opsEnd = run_s2p.run_s2p(ops=ops, db=db)            

def s2p_clean( runops ):

    print('- - - - - - - - -')
    print('S2P CLEANUP')
    ftemp_tmp = runops.get('ftemp')

    sys.stdout.write("remove " + ftemp_tmp + " ?")
    yes = {'yes','y', 'ye', ''}
    no = {'no','n'}
    choice = input().lower()

    if choice in yes:
        try:
            rmtree(ftemp_tmp)
            print('removed ' + ftemp_tmp)
        except:
            print("folder not found or files in use")
    if runops.get('concat'):
        ftiff = runops.get('ftiff')

        sys.stdout.write("remove " + ftiff + " ?")
        choice = input().lower()

        if choice in yes:
            try:
                rmtree(ftiff)
                print('removed ' + ftiff)
            except:
                print("folder not found or files in use")

def si_info( files_fullpath ):
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

    si_info = {
    'dims': mov_dim,
    'level': level,
    'zoom': zoom,
    'framerate': framerate,
    'channels': channels,
    'volumes': volumes,
    'frames': frames,
    'frames_per_file': frames_per_file,
    }

    return si_info

def s2p_ownops( runops ):
    db = {
          'h5py': [], # a single h5 file path
          'h5py_key': 'data', # file paths
          'look_one_level_down': False, # whether to look in ALL subfolders when searching for tiffs
          'data_path': [runops.get('froot')],  # a list of folders with tiffs
                                 # (or folder of folders with tiffs if look_one_level_down is True, or subfolders is not empty)
          'delete_bin': False,
          'save_path0': runops.get('fsave'),
          'fast_disk': runops.get('ftemp'),
          'subfolders': [], # choose subfolders of 'data_path' to look in (optional)

          # main settings
          'nplanes': runops.get('level'), # each tiff has these many planes in sequence
          'nchannels': runops.get('channels'), # each tiff has these many channels per plane
          'functional_chan': 1, # this channel is used to extract functional ROIs (1-based)
          'tau':  1., # this is the main parameter for deconvolution
          'fs': runops.get('framerate') / runops.get('level'),  # sampling rate (PER PLANE - e.g. if you have 12 planes then this should be around 2.5)
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
          'do_registration': 1, # whether to register data (2 forces re-registration)
          'two_step_registration:': True, #default: False) whether or not to run registration twice (for low SNR data). keep_movie_raw must be True for this to work.
          'keep_movie_raw': True,
          'nimg_init': 400, # subsampled frames for finding reference image
          'batch_size': 800, # number of frames per batch
          'maxregshift': 0.05, # max allowed registration shift, as a fraction of frame max(width and height)
          'align_by_chan' : 1, # when multi-channel, you can align by non-functional channel (1-based)
          'reg_tif': False, # whether to save registered tiffs
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
          'tiff_list': runops.get('files_name') # list of tiffs in folder * data_path *!

      }
    return db

def s2p_out_convert( runops ):
    "(Batch) run to convert s2p outout files to our standard ROI / ROIdata mat files"

    # populates s2p ops with the default options

    # Unpack option list
    exp = runops.get('exp')
    concat = runops.get('concat')
    main_root = runops.get('main_root')

    adata = runops.get('adata')
    ftiff = runops.get('ftiff')
    ftemp = runops.get('ftemp')

    darkframes = runops.get('darkframes')
    readfiles = runops.get('readfiles')

    # make subfolders
    ftemp_tmp = os.path.join(ftemp, 'suite2p')
    Path(ftemp_tmp).mkdir(parents=True, exist_ok=True)
    print('- - - - - - - -')
    print('Made directory (ftemp_tmp) ' + ftemp_tmp)
