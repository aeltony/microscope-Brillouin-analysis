#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 7 17:15:35 2021

@author: ame38
"""
import h5py
import os
import numpy as np
import DataFitting

#### Import all data from hdf5 file (output by DAQ software)
# rawData is dictionary of raw data and parameters used
# procData is dictionary of processed data (i.e. fitted + calibrated spectra)
def importData(filename, verbose=False):
    # Create text file for saving exp notes
    name = os.path.splitext(os.path.basename(filename))[0]
    path = os.path.dirname(filename) + '/'
    f_notes = open(path + name + '/' + name + '.txt', 'w')
    f_notes.write('Data Summary: \n')
    # Open data file and load hdf5 structure to dictionary
    f = h5py.File(filename, 'r')
    rawData = {}
    procData = {}
    win_fac = 0 # Number of points either side of peak to include in fit, 0 = no window
    for e in list(f.keys()):
        rawData[e] = {}
        procData[e] = {}
        if f[e].attrs.__contains__('name'):
            rawData[e]['name'] = f[e].attrs['name']
            print(e + ': ' + f[e].attrs['name'])
            f_notes.write('\n' + e + ': ' + f[e].attrs['name'] + '\n')
        scans = list(f[e].keys())
        scanNums = [s[5:] for s in scans]
        scanNums.sort(key = int)
        for n in scanNums:
            s = 'Scan_' + str(n)
            rawData[e][s] = {}
            procData[e][s] = {}
            rawData[e][s]['attrs'] = dict(f[e][s].attrs.items())
            for k in list(f[e][s].keys()):
                rawData[e][s][k] = np.array(f[e][s][k])
            if rawData[e][s]['attrs']['note']:
                print('Processing ' + e + '/' + s + ', ' + rawData[e][s]['attrs']['note'])
                f_notes.write(rawData[e][s]['attrs']['timestamp'] + ' - ' + \
                              s + ': ' + rawData[e][s]['attrs']['note'] + '\n')
            else:
                print('Processing ' + e + '/' + s)
                f_notes.write(rawData[e][s]['attrs']['timestamp'] + ' - ' + s + '\n')
            # Separate sample vs. calibration frames
            calFrames = rawData[e][s]['CalFreq'].shape[1]
            frames = np.array([rawData[e][s]['attrs']['paramtree/Frame Number/X'],
                               rawData[e][s]['attrs']['paramtree/Frame Number/Y'],
                               rawData[e][s]['attrs']['paramtree/Frame Number/Z']])
            procData[e][s]['RawSpecList'] = np.copy(rawData[e][s]['SpecList'])
            procData[e][s]['CalSpecList'] = np.copy(procData[e][s]['RawSpecList'])
            for i in range(calFrames,0,-1):
                procData[e][s]['RawSpecList'] = np.delete(procData[e][s]['RawSpecList'], \
                                                          np.s_[frames[0]::frames[0]+i], 0)             
            for i in range(frames[0],0,-1):
                procData[e][s]['CalSpecList'] = np.delete(procData[e][s]['CalSpecList'], \
                                                          np.s_[::i+calFrames], 0)
            # Fitting Brillouin spectra
            procData[e][s]['IntegrPhotonsList'] = np.zeros(procData[e][s]['RawSpecList'].shape[0])
            procData[e][s]['FreqList'] = np.zeros(procData[e][s]['RawSpecList'].shape[0])
            procData[e][s]['LinewidthList'] = np.zeros(procData[e][s]['RawSpecList'].shape[0])
            procData[e][s]['FreqList_sig'] = np.zeros(procData[e][s]['RawSpecList'].shape[0])
            procData[e][s]['LinewidthList_sig'] = np.zeros(procData[e][s]['RawSpecList'].shape[0])
            procData[e][s]['FittedSpect'] = np.empty(procData[e][s]['RawSpecList'].shape)
            procData[e][s]['FittedCalSpect'] = np.empty(procData[e][s]['CalSpecList'].shape)
            # Find SD / FSR for every (y, z) coordinate
            procData[e][s]['PxDist'] = np.empty(rawData[e][s]['CalFreq'].shape)
            procData[e][s]['PxDist_sig'] = np.empty(rawData[e][s]['CalFreq'].shape)
            procData[e][s]['SDcal'] = np.empty([frames[1]*frames[2]])
            procData[e][s]['FSRcal'] = np.empty([frames[1]*frames[2]])
            procData[e][s]['SDcal_sig'] = np.empty([frames[1]*frames[2]])
            procData[e][s]['FSRcal_sig'] = np.empty([frames[1]*frames[2]])
            for i in range(frames[1]*frames[2]):
                for j in range(calFrames):
                    procData[e][s]['PxDist'][i, j], width, \
                        procData[e][s]['PxDist_sig'][i, j], width_sig, \
                        procData[e][s]['FittedCalSpect'][i*calFrames+j] = \
                        DataFitting.fitSpectrum(np.copy(procData[e][s]['CalSpecList'][i*calFrames+j]),
                                                win_fac, 0.9, 1e-7, 1e-7, verbose)  # 0.9
                try:
                    procData[e][s]['SDcal'][i], procData[e][s]['FSRcal'][i], \
                        procData[e][s]['SDcal_sig'][i], procData[e][s]['FSRcal_sig'][i] = \
                        DataFitting.fitCalCurve(np.copy(procData[e][s]['PxDist'][i]), \
                                                np.copy(rawData[e][s]['CalFreq'][i]), \
                                                1e-7, 1e-7, verbose)
                    # print('[importData] Fitted SD (line %d) = %.3f ± %.5f GHz/px' \
                    #       %(i, procData[e][s]['SDcal'][i], procData[e][s]['SDcal_sig'][i]))
                    # print('[importData] Fitted FSR (line %d) = %.3f ± %.3f GHz' \
                    #       %(i, procData[e][s]['FSRcal'][i], procData[e][s]['FSRcal_sig'][i]))
                except:
                    procData[e][s]['SDcal'][i] = np.nan
                    procData[e][s]['FSRcal'][i] = np.nan
                    if verbose:
                        print('[importData] Could not fit SD + FSR (line %d)' % i)
                # Use SD + FSR to calculate Brillouin frequency shifts
                for j in range(frames[0]):
                    sline = np.copy(procData[e][s]['RawSpecList'][i*frames[0]+j])
                    sline = np.transpose(sline)
                    procData[e][s]['IntegrPhotonsList'][i*frames[0]+j] = \
                        (0.45/0.5)*np.sum(sline)
                    interPeakDist, width, interPeakDist_sig, width_sig, \
                        procData[e][s]['FittedSpect'][i*frames[0]+j] = \
                        DataFitting.fitSpectrum(sline, win_fac, 0.01, 1e-7, 1e-7, verbose)  # 0.1
                    procData[e][s]['FreqList'][i*frames[0]+j] = \
                        0.5*(procData[e][s]['FSRcal'][i] - \
                        procData[e][s]['SDcal'][i]*interPeakDist)
                    procData[e][s]['FreqList_sig'][i*frames[0]+j] = \
                        0.5*np.sqrt( procData[e][s]['FSRcal_sig'][i]**2 + \
                        (interPeakDist*procData[e][s]['SDcal_sig'][i])**2 + \
                        (interPeakDist_sig*procData[e][s]['SDcal'][i])**2)
                    procData[e][s]['LinewidthList'][i*frames[0]+j] = \
                        procData[e][s]['SDcal'][i]*width
                    procData[e][s]['LinewidthList_sig'][i*frames[0]+j] = \
                        np.sqrt((procData[e][s]['SDcal'][i]*width_sig)**2 +\
                        (procData[e][s]['SDcal_sig'][i]*width)**2)

            # Create XYZ array
            procData[e][s]['IntegrPhotonsArr'] = np.reshape(procData[e][s]['IntegrPhotonsList'], (frames[2], frames[1], frames[0]))
            procData[e][s]['FreqArr'] = np.reshape(procData[e][s]['FreqList'], (frames[2], frames[1], frames[0]))
            procData[e][s]['FreqArr_sig'] = np.reshape(procData[e][s]['FreqList_sig'], (frames[2], frames[1], frames[0]))
            procData[e][s]['LWArr'] = np.reshape(procData[e][s]['LinewidthList'], (frames[2], frames[1], frames[0]))
            procData[e][s]['LWArr_sig'] = np.reshape(procData[e][s]['LinewidthList_sig'], (frames[2], frames[1], frames[0]))
            # Match S-shaped line scan pattern
            # procData[e][s]['FreqArr'][1::2, :] = \
            #     np.fliplr(procData[e][s]['FreqArr'][1::2, :])
    f_notes.close()
    f.close()
    return (rawData, procData)