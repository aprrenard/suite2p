import os
import sys
import string

import numpy as np
import PySimpleGUI as sg
import scipy.io as sio

from suite2p import run_s2p

# Populates ops with the default options.
ops = run_s2p.default_ops()
db = []

while True:
    # Solve folders.
    base_folder = '/mnt/synobrice01/Data/2P_Imaging/ANALYSIS'
    event, values = sg.Window('Select folder in analysis directory').Layout(
        [[sg.In(), sg.FolderBrowse(initial_folder=base_folder)],
        [sg.CloseButton('Select'), sg.CloseButton('Cancel')]]
    ).Read()
    analysis_folder = values[0]

    # Load experiment info.
    hinfo_file = os.path.join(analysis_folder, 'ExperimentInfo_suite2p.mat')
    if not os.path.exists(hinfo_file):
        raise FileNotFoundError(
            'No ExperimentInfo_suite2p.mat in {}.\nRun suite2p_aligntrials.m first.'
            .format(analysis_folder))
    hinfo = sio.loadmat(hinfo_file)
    mesc_files = [hinfo['fmescfull'][0][i][0] for i in range(hinfo['fmescfull'][0].shape[0])]
    mesc_fidx = np.array([hinfo['hmesc'][0][i][0][0][0] - 1
                          for i in range(hinfo['hmesc'][0].shape[0])])
    keys = [hinfo['keys'][0][i][0] for i in range(hinfo['keys'][0].shape[0])]

    # Solve compartment.
    event, values = sg.Window('Somas or axons?').Layout(
        [[sg.Text('Were you imaging somas or axons?')],
         [sg.CloseButton('Somas'), sg.CloseButton('Axons')]]
    ).Read()
    if event == 'Somas':
        connected = True
        thr_scaling = 4.0
    elif event == 'Axons':
        connected = False
        thr_scaling = 2.0

    # Solve red channel.
    event, values = sg.Window('Red channel?').Layout(
        [[sg.Text('Was there a red channel for this recording that you wish to process?')],
         [sg.CloseButton('Yes'), sg.CloseButton('No')]]
    ).Read()
    if event == 'Yes':
        nchannels = 2
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

    # Set recording specific parameters.
    idb = {
          'h5py': mesc_files[0], # Analysis folder.
          'look_one_level_down': False, # for h5 files, whether to use all files in same folder
          'data_path': [], # keep this empty!
          'save_folder': analysis_folder,
          'fast_disk': '/home/master',
          'mesc_files': mesc_files,
          'mesc_fidx': mesc_fidx,
          'keys': keys,
          'nchannels': nchannels,
          'reg_tif_chan2': reg_tif_chan2,
          'connected': connected,
          'threshold_scaling': thr_scaling,
          }

    db.append(idb)

    # Batch more recordings?
    event, values = sg.Window('Queue more recordings?').Layout(
        [[sg.Text('Do you want to queue another recording?')],
         [sg.CloseButton('Yes'), sg.CloseButton('No')]]
    ).Read()
    if event == 'No':
        break


# Run suite2p.
error = []
for idb in db:
    try:
        line = ('PROCESSING BATCHED RECORDING {}/{}: {}'
                .format(db.index(idb)+1, len(db), os.path.split(analysis_folder)[1]))
        print(line)
        print('-' * len(line)+'\n')

        opsEnd=run_s2p.run_s2p(ops=ops,db=idb)
        print('\n\n')
    except Exception as e:
        error.append(os.path.split(analysis_folder)[1])
        print(e)
        print("AN ERROR OCCURED. PROCESSING DIDN'T TERMINATE FOR {}.\n"
              .format(os.path.split(analysis_folder)[1]))
        print('\n\n')

# Print recording for which suite2p didn't terminate.
if error:
    print('')
    print('An error occured in {}/{} recordings:\n'.format(len(error), len(db)))
    for e in error:
        print(e)
    print('Check what is wrong and run the script again for those.')
    print('You can open the others in the gui.')
else:
    print('')
    print('Processing successfull for all batched recording!')
    print('You can open them in the gui.')
