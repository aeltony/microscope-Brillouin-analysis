#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 5 19:07:03 2021

@author: ame38
"""
import DataImport
import Plotting
import os

# path = '/Users/ame38/Downloads/2021-06-28/'
path = '/Users/ame38/Dropbox (Partners HealthCare)/Bio-Optics Lab/Data/' + \
        'Microbeads/2022-01-08/'
# path = '/Users/ame38/Dropbox (Partners HealthCare)/Bio-Optics Lab/Data/' + \
#         'Mouse ear/2021-12-20/'
# path = '/Users/ame38/Dropbox (Partners HealthCare)/Bio-Optics Lab/Projects/' + \
#         'Joslin beta cells/2022-03-24/'
# session = 'Control_30hr'
# session = 'Elastase_100ng_28hr'
# session = 'Elastase_250ng_25hr'
# session = 'Elastase_250ng_telaprevir_24hr'
# session = 'Telaprevir_27hr'
# session = 'Whole_mouse_blood'
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
###### Plot images for single exp, single scan ######

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
###### Plot images of all scans for exp ######

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

#%%
###### Process cell group data + plot ######
import CellGroup

# Frequency + linewidth ranges for plotting
fmin = 5.00
fmax = 5.60
lmin = 0.30
lmax = 0.65

# Threshold for masking cell
threshFixed = 0 # 5.11
### If threshFixed = 0, threshold set by Li's iterative Minimum Cross Entropy method

CellGroup.processCellGroup(fullPath, rawData, procData, threshFixed, \
                          fmin, fmax, lmin, lmax)

#%%
###### Process cell group data + plot ######
import CellGroup

# fullPath = '/Users/ame38/Dropbox (Partners HealthCare)/Bio-Optics Lab/Projects/' \
#             + 'Joslin beta cells/2022-02-16/All Groups/'
fullPath = '/Users/ame38/Dropbox (Partners HealthCare)/Bio-Optics Lab/Projects/' \
            + 'Joslin beta cells/2022-04-07/All Groups/'

numCond = 3
numCells = 6
CellGroup.plotAllGroups(fullPath, numCond, numCells)

#%%
###### Process data for all cell dishes + plot ######
import CellGroup

fullPath = '/Users/ame38/Dropbox (Partners HealthCare)/Bio-Optics Lab/Projects/' \
            + 'Joslin beta cells/Aggregate/'

numCond = 5
numCells = 24
CellGroup.plotAllDishes(fullPath, numCond, numCells)

#%%
###### Process before-after data for single cell ######
import CellBeforeAfter

# Frequency + linewidth ranges for plotting
fmin = 5.0
fmax = 5.55
lmin = 0.3
lmax = 0.7

# Cell info
cellNum = 3
exp_1 = 'Exp_0'
scan_1 = 'Scan_2'
exp_2 = 'Exp_2'
scan_2 = 'Scan_2'
threshFixed = 0 # 5.11, 5.22
### If threshFixed = 0, threshold set by Li's iterative Minimum Cross Entropy method

CellBeforeAfter.processBeforeAfter(fullPath, rawData, procData, cellNum, \
                               exp_1, scan_1, exp_2, scan_2, threshFixed, \
                                   fmin, fmax, lmin, lmax)

#%%
###### Process before-after data from txt file ######
import CellBeforeAfter

CellBeforeAfter.plotBeforeAfter(fullPath)

#%%
###### Process aggregate before-after data from txt file ######
import CellBeforeAfter
fullPath = '/Users/ame38/Dropbox (Partners HealthCare)/Bio-Optics Lab/Projects/Joslin beta cells/Aggregate/' + \
                'Telaprevir_before_after.txt'

CellBeforeAfter.plotAggregate(fullPath)

#%%
###### Process aggregate before-after data from txt file ######
import CellBeforeAfter
fullPath = '/Users/ame38/Dropbox (Partners HealthCare)/Bio-Optics Lab/Projects/Joslin beta cells/Aggregate/'

# CellBeforeAfter.plotAllGroups(fullPath)
CellBeforeAfter.plotAllGroupsCorr(fullPath)  # Subtracts medium value for temp. correction

#%%
###### Plot time course for single exp ######
import CellTimeCourse

# Frequency + linewidth ranges for plotting
fmin = 5.0
fmax = 5.55
lmin = 0.3
lmax = 0.7

# Exp to plot
exp = 'Exp_1'
# Threshold for masking cell
threshFixed = 0 # 5.11
### If threshFixed = 0, threshold set by Li's iterative Minimum Cross Entropy method

CellTimeCourse.plotTimeCourse(fullPath, rawData, procData, exp, threshFixed, \
                          fmin, fmax, lmin, lmax)

#%%
import CorneaCrossSection

exp = 'Exp_0'
scan = 'Scan_0'

CorneaCrossSection.processCrossSection(fullPath, rawData, procData, exp, scan)

#%%
import CorneaCrossSection

fullPath = '/Users/ame38/Dropbox (Partners HealthCare)/Bio-Optics Lab/Projects/' + \
                'Anisotropy of cornea/2021 New Data/Cross sections/' + \
                'Porcine_cornea_cross_sections.txt'

CorneaCrossSection.plotAggregate(fullPath)
