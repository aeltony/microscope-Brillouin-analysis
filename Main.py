#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 5 19:07:03 2021

@author: ame38
"""
import DataImport
import Plotting
import os

path = '/Users/ame38/Downloads/2022-01-08/'
session = 'Microbeads'

# Create folder for session
fullPath = path + session + '/'
if not os.path.exists(fullPath):
    os.makedirs(fullPath)

# Import data and fit spectra
filename = path + session + '.hdf5'
(rawData, procData) = DataImport.importData(filename, verbose=True)

# Create folders to save plots
for exp in list(rawData.keys()):
    if not os.path.exists(fullPath + exp + '/'):
        os.makedirs(fullPath + exp + '/')

#%%
###### Plot XY images for single exp, single scan ######

# Save figures?
saveFig = False

# Frequency, linewidth, signal ranges for plotting
fmin = 4.90    # GHz  #5.00 #9.0 #5.3  #5.58
fmax = 6.60    # GHz  #5.60 #9.6 #6.5 #5.45 #5.82
lmin = 0.20    # GHz  #0.30 #0.2
lmax = 2.00    # GHz  #0.65 #0.7
pmin = 1.4e4   # Brillouin photons (integrated) #2.2e4
pmax = 1.1e5   # Brillouin photons (integrated) #3.2e4

# Exp to plot
exp = 'Exp_0'
scan = 'Scan_1'
surfZ = 0 # Z-coordinate of surface

# Plot Brillouin frequency shift, linewidth, and brightfield image
Plotting.plotImages(fullPath, rawData, procData, exp, scan, surfZ, \
                    fmin, fmax, lmin, lmax, pmin, pmax, saveFig)

#%%
###### Plot XY images of all scans for exp ######

# Save figures?
saveFig = True

# Frequency, linewidth, signal ranges for plotting
fmin = 5.00    # GHz  #5.00 #9.0 #5.3  #5.58
fmax = 5.60    # GHz  #5.60 #9.6 #6.5 #5.45 #5.82
lmin = 0.30    # GHz  #0.30 #0.2
lmax = 0.65    # GHz  #0.65 #0.7
pmin = 2.2e4   # Brillouin photons (integrated) #2.2e4
pmax = 3.2e4   # Brillouin photons (integrated) #3.2e4

# Exp to plot
exp = 'Exp_0'
scans = list(rawData[exp].keys())
if scans[0]=='name':
    scans = scans[1:]
surfZ = 0 # Z-coordinate of surface

for scan in scans:
    # Plot Brillouin frequency shift, linewidth, and brightfield image
    Plotting.plotImages(fullPath, rawData, procData, exp, scan, surfZ, \
                        fmin, fmax, lmin, lmax, pmin, pmax, saveFig)

#%%
###### Plot XZ image for single exp, single scan ######

# Save figures?
saveFig = False

# Frequency, linewidth, signal ranges for plotting
fmin = 4.95    # GHz  #5.00 #9.0 #5.3  #5.58
fmax = 5.45    # GHz  #5.60 #9.6 #6.5 #5.45 #5.82
lmin = 0.20    # GHz  #0.30 #0.2
lmax = 0.70    # GHz  #0.65 #0.7
pmin = 1.0e4   # Brillouin photons (integrated) #2.2e4
pmax = 0.8e5   # Brillouin photons (integrated) #3.2e4

# Exp to plot
exp = 'Exp_0'
scan = 'Scan_0'
surfZ = 0 # Z-coordinate of surface

# Plot Brillouin frequency shift, linewidth, and brightfield image
Plotting.plotXZimage(fullPath, rawData, procData, exp, scan, \
                    fmin, fmax, lmin, lmax, pmin, pmax, saveFig)

#%%
###### Plot Brightfield images (only) of all scans for exp ######

# Save figures?
saveFig = False

# Exp to plot
exp = 'Exp_11'
scans = list(rawData[exp].keys())
if scans[0]=='name':
    scans = scans[1:]
surfZ = 0 # Z-coordinate of surface

for scan in scans:
    # Plot Brillouin frequency shift, linewidth, and brightfield image
    Plotting.plotBrightfield(fullPath, rawData, procData, exp, scan, surfZ, saveFig)

#%%
###### Plot Brillouin spectra (data) for a single line of a scan ######
exp = 'Exp_3'
scan = 'Scan_9'
frame = 0  # Frame is the z-stack number
line = 0  # Lines correspond to fixed y positions

# Save figures?
saveFig = False
if saveFig:
    # Create folder to save plots:
    if not os.path.exists(fullPath + exp + '/' + scan + '/'):
        os.makedirs(fullPath + exp + '/' + scan + '/')

Plotting.plotSpectra(fullPath, rawData, procData, exp, scan, frame, line, saveFig)

# # Plot all frames at once
# for f in range(0,200):
#     frame = f
#     Plotting.plotSpectra(fullPath, rawData, procData, exp, scan, frame, line, saveFig)

#%%
###### Plot raw spectrum image (from Andor camera) for a single line of a scan ######
exp = 'Exp_3'
scan = 'Scan_9'
frame = 0  # Frame is the z-stack number
line = 0  # Lines correspond to fixed y positions

# Save figures?
saveFig = False
if saveFig:
    # Create folder to save plots:
    if not os.path.exists(fullPath + exp + '/' + scan + '/'):
        os.makedirs(fullPath + exp + '/' + scan + '/')

Plotting.plotAndorImage(fullPath, rawData, exp, scan, frame, line, saveFig)

#%%
###### Plot calibration spectra for a single line of a scan ######
exp = 'Exp_0'
scan = 'Scan_0'
frame = 0 # Frame is the z-stack number
line = 0 # Lines correspond to fixed y positions

# Save figures?
saveFig = False
if saveFig:
    # Create folder to save plots:
    if not os.path.exists(fullPath + exp + '/' + scan + '/'):
        os.makedirs(fullPath + exp + '/' + scan + '/')

Plotting.plotCalSpectra(fullPath, rawData, procData, exp, scan, frame, line, saveFig)

# # Plot all frames at once
# for f in range(0,200):
#     frame = f
#     Plotting.plotCalSpectra(fullPath, rawData, procData, exp, scan, frame, line, saveFig)

#%%
###### Plot calibration spectrum image (from Andor camera) for a single line of a scan ######
exp = 'Exp_0'
scan = 'Scan_0'
frame = 0  # Frame is the z-stack number
line = 0  # Lines correspond to fixed y positions

# Save figures?
saveFig = False
if saveFig: 
    # Create folder to save plots:
    if not os.path.exists(fullPath + exp + '/' + scan + '/'):
        os.makedirs(fullPath + exp + '/' + scan + '/')

Plotting.plotCalAndorImage(fullPath, rawData, exp, scan, frame, line, saveFig)
