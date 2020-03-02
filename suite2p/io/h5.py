import gc
import numpy as np
from natsort import natsorted
import math, time
import glob, h5py, os, json
from scipy import signal
from scipy import stats, signal
from scipy.sparse import linalg
import scipy.io
from . import utils

def h5py_to_binary(ops):
    """  finds h5 files and writes them to binaries

    Parameters
    ----------
    ops : dictionary
        'nplanes', 'h5_path', 'h5_key', 'save_path', 'save_folder', 'fast_disk',
        'nchannels', 'keep_movie_raw', 'look_one_level_down'

    Returns
    -------
        ops1 : list of dictionaries
            'Ly', 'Lx', ops1[j]['reg_file'] or ops1[j]['raw_file'] is created binary

    """
    ops1 = utils.init_ops(ops)

    nplanes = ops1[0]['nplanes']
    nchannels = ops1[0]['nchannels']

    # open all binary files for writing
    if ops['mesc']:
        ops1, h5list, reg_file, reg_file_chan2 = utils.find_files_open_binaries(ops1, True)
        for ops in ops1:
            if 'data_path' in ops and len(ops['data_path'])==0:
                ops['data_path'] = [os.path.dirname(ops['h5py'])]
            elif 'data_path' not in ops:
                ops['data_path'] = [os.path.dirname(ops['h5py'])]
        ops1[0]['h5list'] = ops1[0]['mesc_files']
        mesc_files = ops1[0]['mesc_files']
        mesc_fidx = ops1[0]['mesc_fidx']
        keys = ops1[0]['keys']

        iall = 0
        for j in range(ops['nplanes']):
            ops1[j]['nframes_per_folder'] = np.zeros(len(h5list), np.int32)
    else:
        ops1, h5list, reg_file, reg_file_chan2 = utils.find_files_open_binaries(ops1, True)
        for ops in ops1:
            if 'data_path' in ops and len(ops['data_path'])==0:
                ops['data_path'] = [os.path.dirname(ops['h5py'])]
            elif 'data_path' not in ops:
                ops['data_path'] = [os.path.dirname(ops['h5py'])]
        ops1[0]['h5list'] = h5list
        keys = ops1[0]['h5py_key']
        if isinstance(keys, str):
            keys = [keys]
        iall = 0
        for j in range(ops['nplanes']):
            ops1[j]['nframes_per_folder'] = np.zeros(len(h5list), np.int32)

    if ops['mesc']:
        for ikey in keys:

            ifileidx = mesc_fidx[keys.index(ikey)]
            ifile = mesc_files[ifileidx]
            with h5py.File(ifile, 'r') as f:
                # Reading from mesc file.
                # keep track of the plane identity of the first frame (channel identity is assumed always 0)
                nbatch = nplanes*nchannels*math.ceil(ops1[0]['batch_size']/(nplanes*nchannels))
                print(keys)
                print(mesc_files)
                nframes_all = f[ikey + '/Channel_1'].shape[0]

                if ops1[0]['nchannels'] > 1:
                    # If two channels, data is read from two different
                    # datasets in mesc file.
                    nframes_all = 2*nframes_all
                nbatch = min(nbatch, nframes_all)
                if nchannels>1:
                    nfunc = ops['functional_chan'] - 1
                else:
                    nfunc = 0
                # loop over all tiffs
                ik = 0
                while 1:
                    if ops1[0]['nchannels'] > 1:
                        irange = np.arange(ik, min(ik+nbatch, nframes_all/2), 1)
                    else:
                        irange = np.arange(ik, min(ik+nbatch, nframes_all), 1)
                    if irange.size==0:
                        break

                    if ops1[0]['nchannels'] > 1:
                        key_1 = ikey + '/Channel_0'
                        key_2 = ikey + '/Channel_1'
                        im_1 = f[key_1][irange, :, :]
                        im_1 = np.invert(im_1)
                        im_2 = f[key_2][irange, :, :]
                        im_2 = np.invert(im_2)
                        # Interleave the channels.
                        s = im_1.shape
                        im = np.zeros((s[0]*2, s[1], s[2]))
                        im[0::2] = im_1
                        im[1::2] = im_2
                    else:
                        im = f[ikey + '/channel_0']
                        im = np.invert(im)

                        nframes = im.shape[0]
                        if type(im[0,0,0]) == np.uint16:
                            im = im / 2
                        for j in range(0,nplanes):
                            if iall==0:
                                ops1[j]['meanImg'] = np.zeros((im.shape[1],im.shape[2]),np.float32)
                                if nchannels>1:
                                    ops1[j]['meanImg_chan2'] = np.zeros((im.shape[1],im.shape[2]),np.float32)
                                ops1[j]['nframes'] = 0
                            i0 = nchannels * ((j)%nplanes)
                            im2write = im[np.arange(int(i0)+nfunc, nframes, nplanes*nchannels),:,:].astype(np.int16)
                            reg_file[j].write(bytearray(im2write))
                            ops1[j]['meanImg'] += im2write.astype(np.float32).sum(axis=0)
                            if nchannels>1:
                                im2write = im[np.arange(int(i0)+1-nfunc, nframes, nplanes*nchannels),:,:].astype(np.int16)
                                reg_file_chan2[j].write(bytearray(im2write))
                                ops1[j]['meanImg_chan2'] += im2write.astype(np.float32).sum(axis=0)
                            ops1[j]['nframes'] += im2write.shape[0]
                            ops1[j]['nframes_per_folder'][ifileidx] += im2write.shape[0]
                        if ops1[0]['nchannels'] > 1:
                            # Because need to take half im.shape[0] frames from
                            # each of the two channels.
                            ik += int(nframes/2)
                            iall += int(nframes/2)
                        else:
                            ik += nframes
                            iall += nframes
    else:
        for ih5,h5 in enumerate(h5list):
            with h5py.File(h5, 'r') as f:
                # if h5py data is 4D instead of 3D, assume that
                # data = nframes x nplanes x pixels x pixels
                for key in keys:
                    # Reading from mesc file.
                    hdims = f[key].ndim
                    # keep track of the plane identity of the first frame (channel identity is assumed always 0)
                    nbatch = nplanes*nchannels*math.ceil(ops1[0]['batch_size']/(nplanes*nchannels))
                    if hdims==3:
                        nframes_all = f[key].shape[0]
                    else:
                        nframes_all = f[key].shape[0] * f[key].shape[1]
                    nbatch = min(nbatch, nframes_all)
                    if nchannels>1:
                        nfunc = ops['functional_chan'] - 1
                    else:
                        nfunc = 0
                    # loop over all tiffs
                    ik = 0
                    while 1:
                        if hdims==3:
                            irange = np.arange(ik, min(ik+nbatch, nframes_all), 1)
                            if irange.size==0:
                                break
                            im = f[key][irange, :, :]
                        else:
                            irange = np.arange(ik/nplanes, min(ik/nplanes+nbatch/nplanes, nframes_all/nplanes), 1)
                            if irange.size==0:
                                break
                            im = f[key][irange,:,:,:]
                            im = np.reshape(im, (im.shape[0]*nplanes,im.shape[2],im.shape[3]))
                        nframes = im.shape[0]
                        if type(im[0,0,0]) == np.uint16:
                            im = im / 2
                        for j in range(0,nplanes):
                            if iall==0:
                                ops1[j]['meanImg'] = np.zeros((im.shape[1],im.shape[2]),np.float32)
                                if nchannels>1:
                                    ops1[j]['meanImg_chan2'] = np.zeros((im.shape[1],im.shape[2]),np.float32)
                                ops1[j]['nframes'] = 0
                            i0 = nchannels * ((j)%nplanes)
                            im2write = im[np.arange(int(i0)+nfunc, nframes, nplanes*nchannels),:,:].astype(np.int16)
                            reg_file[j].write(bytearray(im2write))
                            ops1[j]['meanImg'] += im2write.astype(np.float32).sum(axis=0)
                            if nchannels>1:
                                im2write = im[np.arange(int(i0)+1-nfunc, nframes, nplanes*nchannels),:,:].astype(np.int16)
                                reg_file_chan2[j].write(bytearray(im2write))
                                ops1[j]['meanImg_chan2'] += im2write.astype(np.float32).sum(axis=0)
                            ops1[j]['nframes'] += im2write.shape[0]
                            ops1[j]['nframes_per_folder'][ih5] += im2write.shape[0]
                        ik += nframes
                        iall += nframes

    # write ops files
    do_registration = ops1[0]['do_registration']
    do_nonrigid = ops1[0]['nonrigid']
    for ops in ops1:
        ops['Ly'] = im2write.shape[1]
        ops['Lx'] = im2write.shape[2]
        if not do_registration:
            ops['yrange'] = np.array([0,ops['Ly']])
            ops['xrange'] = np.array([0,ops['Lx']])
        ops['meanImg'] /= ops['nframes']
        if nchannels>1:
            ops['meanImg_chan2'] /= ops['nframes']
        np.save(ops['ops_path'], ops)
    # close all binary files and write ops files
    for j in range(0,nplanes):
        reg_file[j].close()
        if nchannels>1:
            reg_file_chan2[j].close()
    return ops1
