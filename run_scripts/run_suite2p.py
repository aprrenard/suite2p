import os
import sys
import string

import numpy as np
import PySimpleGUI as sg
import scipy.io as sio

from suite2p import run_s2p

# set your options for running
ops = run_s2p.default_ops() # populates ops with the default options

# Solve folders.
# --------------

base_folder = '/mnt/synobrice01/Data/2P_Imaging/ANALYSIS'
event, values = sg.Window('Select folder in analysis directory').Layout(
    [[sg.In(), sg.FolderBrowse(initial_folder=base_folder)],
    [sg.CloseButton('Select'), sg.CloseButton('Cancel')]]
).Read()
analysis_folder = values[0]


# Solve red channel.
# ------------------

event, values = sg.Window('Red channel?').Layout(
    [[sg.Text('Was there a red channel for this recording (channel 2)?')],
     [sg.CloseButton('Yes'), sg.CloseButton('No')]]
).Read()

if event == 'Yes':
    nchannels = 2

    event, values = sg.Window('Registration on red channel?').Layout(
        [[sg.Text('Do you want to do the movement correction on the red channel?')],
         [sg.CloseButton('Yes'), sg.CloseButton('No')]]
    ).Read()
    if event == 'Yes':
        align_by_chan = True
    else:
        align_by_chan = False

    event, values = sg.Window('Save registered red tiff?').Layout(
        [[sg.Text('Do you want to create movement corrected tiff files for the red channel?')],
         [sg.CloseButton('Yes'), sg.CloseButton('No')]]
    ).Read()
    if event == 'Yes':
        reg_tif_chan2 = True
    else:
        reg_tif_chan2 = False

else:
    nchannels = 1
    reg_tif_chan2 = False
    align_by_chan = False


# Load experiment info.
# ---------------------

hinfo_file = os.path.join(analysis_folder, 'ExperimentInfo_suite2p.mat')
hinfo = sio.loadmat(hinfo_file)
mesc_files = [hinfo['fmescfull'][0][i][0] for i in range(hinfo['fmescfull'][0].shape[0])]
mesc_fidx = np.array([hinfo['hmesc'][0][i][0][0][0] - 1
                      for i in range(hinfo['hmesc'][0].shape[0])])
keys = [hinfo['keys'][0][i][0] for i in range(hinfo['keys'][0].shape[0])]


# Set parameters.
# ---------------

db = {
      'h5py': mesc_files[0], # Analysis folder.
      'look_one_level_down': False, # for h5 files, whether to use all files in same folder
      'data_path': [], # keep this empty!
      'save_folder': analysis_folder,
      'fast_disk': '/home/master',
      'mesc_files': mesc_files,
      'mesc_fidx': mesc_fidx,
      'keys': keys,
      'nchannels': nchannels,
      'red_tif_chan2': reg_tif_chan2,
      'align_by_chan': align_by_chan,

      'delete_bin': False,
      'sparse_mode': True,
      'spatial_scale': 1,
      'threshold_scaling': 4.0
    }

# run one experiment
opsEnd=run_s2p.run_s2p(ops=ops,db=db)
