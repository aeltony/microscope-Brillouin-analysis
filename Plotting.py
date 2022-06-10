#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 23 11:07:33 2021

@author: ame38
"""
import os
import matplotlib.pyplot as plt
plt.rcParams['figure.max_open_warning'] = False
plt.rcParams['pdf.fonttype'] = 42
from matplotlib.colors import LinearSegmentedColormap
from matplotlib_scalebar.scalebar import ScaleBar
import numpy as np
plt.style.use('myPlotStyle.mplstyle')


#### Plot raw spectra images from Andor camera
def plotAndorImage(path, rawData, exp, scan, frame, line, saveFig=False):
    if saveFig:
        # Make folder to save images
        if not os.path.exists(path + exp + '/' + scan + '/' + 'Spectra' + '/'):
            os.makedirs(path + exp + '/' + scan + '/'  + 'Spectra' + '/')
    steps = rawData[exp][scan]['attrs']['paramtree/Frame Number/X']
    lines = rawData[exp][scan]['attrs']['paramtree/Frame Number/Y']
    calFrames = rawData[exp][scan]['CalFreq'].shape[1]
    for s in range(steps):
        idx = frame*lines*(steps + calFrames) + line*(steps + calFrames) + s
        image = rawData[exp][scan]['AndorImage'][idx]
        plt.figure()
        plt.grid(b=None)
        plt.title(exp + '/' + scan + ': frame %d, line %d, step %d'\
                  %(frame, line, s))
        plt.imshow(image, cmap='Greys_r', vmin=np.amin(image), vmax=np.amax(image), \
                    interpolation='nearest')
        ax = plt.gca()
        ax.axis('off')
        ax.autoscale_view

        if saveFig:
            plt.savefig(path + exp + '/' + scan + '/' + 'Spectra' + '/' \
                + exp + '_' + scan + r'_frame_%d_line_%d_step_%d_raw.png' \
                %(frame, line, s), \
                transparent=True, bbox_inches='tight', pad_inches=0.01)


#### Plot calibration spectra images from Andor camera
def plotCalAndorImage(path, rawData, exp, scan, frame, line, saveFig=False):
    if saveFig:
        # Make folder to save images
        if not os.path.exists(path + exp + '/' + scan + '/' + 'CalSpectra' + '/'):
            os.makedirs(path + exp + '/' + scan + '/'  + 'CalSpectra' + '/')
    steps = rawData[exp][scan]['attrs']['paramtree/Frame Number/X']
    lines = rawData[exp][scan]['attrs']['paramtree/Frame Number/Y']
    calFrames = rawData[exp][scan]['CalFreq'].shape[1]
    for c in range(calFrames):
        idx = frame*lines*(steps + calFrames) + line*(steps + calFrames) + steps + c
        image = rawData[exp][scan]['AndorImage'][idx]
        plt.figure()
        plt.grid(b=None)
        plt.title(exp + '/' + scan + ': frame %d, line %d, calFrame %d'\
                  %(frame, line, c))
        plt.imshow(image, cmap='Greys_r', vmin=np.amin(image), vmax=np.amax(image), \
                    interpolation='nearest')
        ax = plt.gca()
        ax.axis('off')
        ax.autoscale_view

        if saveFig:
            plt.savefig(path + exp + '/' + scan + '/' + 'CalSpectra' + '/' \
                + exp + '_' + scan + r'_frame_%d_line_%d_calFrame_%d_raw.png' \
                %(frame, line, c), \
                transparent=True, bbox_inches='tight', pad_inches=0.01)


#### Plot Brillouin spectra (data)
def plotSpectra(path, rawData, procData, exp, scan, frame, line, saveFig=False):
    if saveFig:
        # Make folder to save images
        if not os.path.exists(path + exp + '/' + scan + '/' + 'Spectra' + '/'):
            os.makedirs(path + exp + '/' + scan + '/'  + 'Spectra' + '/')
    steps = rawData[exp][scan]['attrs']['paramtree/Frame Number/X']
    lines = rawData[exp][scan]['attrs']['paramtree/Frame Number/Y']
    for s in range(steps):
        idx = frame*lines*steps + line*steps + s
        plt.figure()
        plt.plot(procData[exp][scan]['FittedSpect'][idx], \
                 '-', linewidth=2)
        plt.plot(procData[exp][scan]['RawSpecList'][idx], \
                 'ko', markersize=4)
        plt.title(exp + '/' + scan + ': frame %d, line %d, step %d \n Fitted BFS = %.3f ± %.3f GHz \n Fitted LW = %.3f ± %.3f GHz'\
                  %(frame, line, s, procData[exp][scan]['FreqList'][idx], \
                  procData[exp][scan]['FreqList_sig'][idx], \
                  procData[exp][scan]['LinewidthList'][idx], \
                  procData[exp][scan]['LinewidthList_sig'][idx]))
        plt.xlabel('Pixels')
        plt.ylabel('Counts')
        if saveFig:
            plt.savefig(path + exp + '/' + scan + '/' + 'Spectra' + '/' \
                + exp + '_' + scan + r'_frame_%d_line_%d_step_%d.png' \
                %(frame, line, s), \
                transparent=True, bbox_inches='tight', pad_inches=0.01)


#### Plot calibration spectra
def plotCalSpectra(path, rawData, procData, exp, scan, frame, line, saveFig=False):
    if saveFig:
        # Make folder to save images
        if not os.path.exists(path + exp + '/' + scan + '/' + 'CalSpectra' + '/'):
            os.makedirs(path + exp + '/' + scan + '/'  + 'CalSpectra' + '/')
    lines = rawData[exp][scan]['attrs']['paramtree/Frame Number/Y']
    calFrames = rawData[exp][scan]['CalFreq'].shape[1]
    for s in range(calFrames):
        plt.figure()
        plt.plot(procData[exp][scan]['FittedCalSpect'][frame*lines*calFrames + line*calFrames + s], \
                 '-', linewidth=2)
        plt.plot(procData[exp][scan]['CalSpecList'][frame*lines*calFrames + line*calFrames + s], \
                 'ko', markersize=4)
        plt.title(exp + '/' + scan + ': frame %d, line %d, calFrame %d \n EOM setting: %.3f GHz' \
            %(frame, line, s, rawData[exp][scan]['CalFreq'][line, s]))
        plt.xlabel('Pixels')
        plt.ylabel('Counts')
        if saveFig:
            plt.savefig(path + exp + '/' + scan + '/' + 'CalSpectra' + '/' \
                + exp + '_' + scan + r'_frame_%d_line_%d_calFrame_%d.png' %(frame, line, s), \
                transparent=True, bbox_inches='tight', pad_inches=0.01)


#### Plot all XY images for single exp and scan
def plotImages(path, rawData, procData, exp, scan, surfZ=0, fmin=5.58, fmax=5.82, lmin=0.0, lmax=1.5, pmin=0.0, pmax=1e5, saveFig=False):  
    # Import standard Brillouin colormap
    colors = np.genfromtxt('colormap.txt', delimiter='\t')
    colormap = LinearSegmentedColormap.from_list('brillouin', colors, N=200)
    # colormap = 'brillouin_r'
    # colormap = 'jet'
    # colormap = 'plasma' # perceptually uniform
    
    if saveFig:        
        # Make folder to save PDFs
        if not os.path.exists(path + exp + '/PDFs/'):
            os.makedirs(path + exp + '/PDFs/')
    
    frames = np.array([rawData[exp][scan]['attrs']['paramtree/Frame Number/X'],
                        rawData[exp][scan]['attrs']['paramtree/Frame Number/Y'],
                        rawData[exp][scan]['attrs']['paramtree/Frame Number/Z']])
    step = np.array([rawData[exp][scan]['attrs']['paramtree/Step Size/X'],
                        rawData[exp][scan]['attrs']['paramtree/Step Size/Y'],
                        rawData[exp][scan]['attrs']['paramtree/Step Size/Z']])

    # Measured step is based on the coordinates recorded by the Zaber stage
    measStep = np.array([0.0, 0.0, 0.0])
    if frames[0]>1:
        measStep[0] = rawData[exp][scan]['MotorCoords'][1, 0] - rawData[exp][scan]['MotorCoords'][0, 0]
    if frames[1]>1:
        measStep[1] = rawData[exp][scan]['MotorCoords'][frames[0], 1] - rawData[exp][scan]['MotorCoords'][0, 1]
    if frames[2]>1:
        measStep[2] = rawData[exp][scan]['MotorCoords'][frames[0]*frames[1], 2] - rawData[exp][scan]['MotorCoords'][0, 2]
    
    # Extent of Brillouin map (based on span of Zaber motor coordinates)
    xrange = rawData[exp][scan]['MotorCoords'][frames[0]-1, 0] - \
        rawData[exp][scan]['MotorCoords'][0, 0]
    yrange = rawData[exp][scan]['MotorCoords'][frames[0]*frames[1]-1, 1] - \
        rawData[exp][scan]['MotorCoords'][0, 1]
    
    # Size of single FLIR camera frame in um
    imageSize = 9000.0/rawData[exp][scan]['attrs']['paramtree/Microscope Camera/Magnification']
    
    # Number of images per depth (X-Y) scan
    if 'CMOSImage' in rawData[exp][scan].keys():
        imPerDepth = int(rawData[exp][scan]['CMOSImage'].shape[0]/frames[2])
    elif 'BrightfieldImage' in rawData[exp][scan].keys():
        imPerDepth = int(rawData[exp][scan]['BrightfieldImage'].shape[0]/frames[2])
    
    for z in range(procData[exp][scan]['FreqArr'].shape[0]):
        if surfZ==0:
            depth = z*measStep[2]
        else:
            depth = z*measStep[2] + rawData[exp][scan]['MotorCoords'][0,2] - surfZ
        print('Frame #%d, depth = %.1f um' %(z, depth))
        
        ### Plot composite Brightfield image:
        if 'CMOSImage' in rawData[exp][scan].keys():
            images = np.rot90(rawData[exp][scan]['CMOSImage'][z*imPerDepth:(z+1)*imPerDepth], -1, (2,1)) # Rotate by -90 deg.
        elif 'BrightfieldImage' in rawData[exp][scan].keys():
            images = np.rot90(rawData[exp][scan]['BrightfieldImage'][z*imPerDepth:(z+1)*imPerDepth], -1, (2,1)) # Rotate by -90 deg.
        pxPerMicron = images.shape[2]/imageSize
        xinterval = int(np.ceil(0.5*imageSize/step[0])*step[0]*pxPerMicron)
        yinterval = int(np.ceil(0.5*imageSize/step[1])*step[1]*pxPerMicron)
        imPerLine = 1 + int(np.floor(frames[0]*step[0]*pxPerMicron/xinterval))
        imPerCol = 1 + int(np.floor(frames[1]*step[1]*pxPerMicron/yinterval))
        image = np.zeros((yinterval*imPerCol + images.shape[2] - yinterval, \
                              xinterval*imPerLine + images.shape[1] - xinterval))
        for line in range(imPerCol):
            for im in range(imPerLine):
                image[yinterval*line:yinterval*line + images.shape[2], \
                      xinterval*im:xinterval*im + images.shape[1]] = \
                    images[line*imPerLine + im]
                    
        plt.figure()
        plt.grid(b=None)
        # plt.xlabel(r'x (${\rm \mu m}$)')
        # plt.ylabel(r'y (${\rm \mu m}$)')
        plt.title(exp + '/' + scan + r', z = %.1f ${\rm \mu m}$' %depth)
        plt.imshow(image, cmap='Greys_r', vmin=np.amin(image), vmax=np.amax(image), \
                    interpolation='nearest',
                    extent=(0, 0.5*imageSize*(imPerLine+1), \
                            0, 0.5*imageSize*(imPerCol+1)))
        ax = plt.gca()
        ax.axis('off')
        ax.autoscale_view
        scalebar = ScaleBar(1, 'um', length_fraction=0.2, location='upper left', \
                            color=(1.0, 0.0, 0.0), box_color='None')
        ax.add_artist(scalebar)
        # ax.xaxis.set_ticks_position('bottom')
        if saveFig:
            plt.savefig(path + exp + '/' + exp + '_' + scan + r'_z_%d.png' %z, \
                    transparent=True, bbox_inches='tight', pad_inches=0.01)
        
        ### Plot zoomed-in Brightfield image
        xStartPx = int(np.floor(0.5*image.shape[0]) + \
                        pxPerMicron*rawData[exp][scan]['attrs']['paramtree/More Settings/Laser Focus X'])
        yStartPx = int(np.floor(0.5*image.shape[1]) - \
                        pxPerMicron*rawData[exp][scan]['attrs']['paramtree/More Settings/Laser Focus Y'])
        # xStartPx = 1093
        # yStartPx = 920
        cropImage = image[yStartPx : yStartPx + int(np.round(1.0*yrange*pxPerMicron)), \
                          xStartPx : xStartPx + int(np.round(1.0*xrange*pxPerMicron))]
        plt.figure()
        plt.grid(b=None)
        plt.imshow(cropImage, cmap='Greys_r', vmin=np.amin(cropImage), \
                   vmax=np.amax(cropImage), interpolation='nearest', \
                   extent=(0, 1.0*xrange, 0, 1.0*yrange))
        ax = plt.gca()
        ax.axis('off')
        ax.autoscale_view
        scalebar = ScaleBar(1, 'um', length_fraction=0.3, location='upper left', \
                            color=(1.0, 0.0, 0.0), box_color='None')
        ax.add_artist(scalebar)
        # ax.xaxis.set_ticks_position('bottom')
        if saveFig:
            plt.savefig(path + exp + '/' + exp + '_' + scan + r'_z_%d_cropped.png' %z, \
                        transparent=True, bbox_inches='tight', pad_inches=0.01)
            plt.savefig(path + exp + '/PDFs/' + exp + '_' + scan + r'_z_%d_cropped.pdf' %z, \
                        transparent=True, bbox_inches='tight', pad_inches=0.01)
        
        ### Plot frequency shift
        plt.matshow(np.flip(procData[exp][scan]['FreqArr'][z],0), cmap=colormap, \
                    vmin=fmin, vmax=fmax, origin='lower', interpolation=None, \
                    extent=(0, xrange, 0, yrange))
        # plt.title(exp + '/' + scan + r', z = %.1f ${\rm \mu m}$' %depth)
        plt.xlabel(r'x (${\rm \mu m}$)')
        plt.ylabel(r'y (${\rm \mu m}$)')
        plt.colorbar(label='Brillouin frequency shift [GHz]', shrink=0.75)
        plt.grid(b=None)
        ax = plt.gca()
        ax.axis('off')
        ax.autoscale_view
        scalebar = ScaleBar(1, 'um', length_fraction=0.3, location='upper left', \
                            color='w', box_color='None')
        ax.add_artist(scalebar)
        # ax.xaxis.set_ticks_position('bottom')
        if saveFig:
            plt.savefig(path + exp + '/' + exp + '_' + scan + '_BFS' + r'_z_%d.png' %z, \
                        transparent=True, bbox_inches='tight', pad_inches=0.01)
            plt.savefig(path + exp + '/PDFs/' + exp + '_' + scan + '_BFS' + r'_z_%d.pdf' %z, \
                        transparent=True, bbox_inches='tight', pad_inches=0.01)
    
        ### Plot linewidth
        plt.matshow(np.flip(procData[exp][scan]['LWArr'][z],0), cmap=colormap, \
                    vmin=lmin, vmax=lmax, origin='lower', interpolation=None, \
                    extent=(0, xrange, 0, yrange))
        # plt.title(exp + '/' + scan + r', z = %.1f ${\rm \mu m}$' %depth)
        plt.xlabel(r'x (${\rm \mu m}$)')
        plt.ylabel(r'y (${\rm \mu m}$)')
        plt.colorbar(label='Brillouin peak linewidth [GHz]', shrink=0.75)
        plt.grid(b=None)
        ax = plt.gca()
        ax.axis('off')
        ax.autoscale_view
        scalebar = ScaleBar(1, 'um', length_fraction=0.3, location='upper left', \
                            color='w', box_color='None')
        ax.add_artist(scalebar)
        # ax.xaxis.set_ticks_position('bottom')
        if saveFig:
            plt.savefig(path + exp + '/' + exp + '_' + scan + '_LW' + r'_z_%d.png' %z, \
                        transparent=True, bbox_inches='tight', pad_inches=0.01)
        
        ### Plot signal intensity (integrated Brillouin photons in S/AS peaks)
        plt.matshow(np.flip(procData[exp][scan]['IntegrPhotonsArr'][z],0), cmap=colormap, \
                    vmin=pmin, vmax=pmax, origin='lower', interpolation=None, \
                    extent=(0, xrange, 0, yrange))
        # plt.title(exp + '/' + scan + r', z = %.1f ${\rm \mu m}$' %depth)
        plt.xlabel(r'x (${\rm \mu m}$)')
        plt.ylabel(r'y (${\rm \mu m}$)')
        plt.colorbar(label='Integrated S/AS photons [#]', shrink=0.75)
        plt.grid(b=None)
        ax = plt.gca()
        ax.axis('off')
        ax.autoscale_view
        scalebar = ScaleBar(1, 'um', length_fraction=0.3, location='upper left', \
                            color='w', box_color='None')
        ax.add_artist(scalebar)
        # ax.xaxis.set_ticks_position('bottom')
        if saveFig:
            plt.savefig(path + exp + '/' + exp + '_' + scan + '_photons' + r'_z_%d.png' %z, \
                        transparent=True, bbox_inches='tight', pad_inches=0.01)        


#### Plot XZ image for single exp and scan
def plotXZimage(path, rawData, procData, exp, scan, fmin=5.58, fmax=5.82, lmin=0.0, lmax=1.5, pmin=0.0, pmax=1e5, saveFig=False):
    # Import standard Brillouin colormap
    # colors = np.genfromtxt('colormap.txt', delimiter='\t')
    # colormap = LinearSegmentedColormap.from_list('brillouin', colors, N=200)
    colormap = 'jet'
    # colormap = 'plasma' # perceptually uniform
    
    if saveFig:
        # Make folder to save PDFs
        if not os.path.exists(path + exp + '/PDFs/'):
            os.makedirs(path + exp + '/PDFs/')
        # Make folder to save brightfield images
        if not os.path.exists(path + exp + '/Brightfield/'):
            os.makedirs(path + exp + '/Brightfield/')
    
    frames = np.array([rawData[exp][scan]['attrs']['paramtree/Frame Number/X'],
                        rawData[exp][scan]['attrs']['paramtree/Frame Number/Y'],
                        rawData[exp][scan]['attrs']['paramtree/Frame Number/Z']])

    # Measured step is based on the coordinates recorded by the Zaber stage
    measStep = np.array([0.0, 0.0, 0.0])
    if frames[0]>1:
        measStep[0] = rawData[exp][scan]['MotorCoords'][1, 0] - rawData[exp][scan]['MotorCoords'][0, 0]
    if frames[1]>1:
        measStep[1] = rawData[exp][scan]['MotorCoords'][frames[0], 1] - rawData[exp][scan]['MotorCoords'][0, 1]
    if frames[2]>1:
        measStep[2] = rawData[exp][scan]['MotorCoords'][frames[0]*frames[1], 2] - rawData[exp][scan]['MotorCoords'][0, 2]
    
    # Extent of Brillouin map (based on span of Zaber motor coordinates)
    xrange = rawData[exp][scan]['MotorCoords'][frames[0]-1, 0] - \
        rawData[exp][scan]['MotorCoords'][0, 0]
    zrange = rawData[exp][scan]['MotorCoords'][frames[0]*frames[2]-1, 2] - \
        rawData[exp][scan]['MotorCoords'][0, 2]
        
    ### Plot frequency shift
    plt.matshow(procData[exp][scan]['FreqArr'][:,0,:], cmap=colormap, \
                vmin=fmin, vmax=fmax, origin='lower', interpolation=None, \
                extent=(0, xrange, 0, zrange))
    plt.colorbar(label='Brillouin frequency shift [GHz]', shrink=0.75)
    plt.grid(b=None)
    ax = plt.gca()
    ax.axis('off')
    ax.autoscale_view
    scalebar = ScaleBar(1, 'um', length_fraction=0.3, location='upper left', \
                        color='w', box_color='None')
    ax.add_artist(scalebar)
    if saveFig:
        plt.savefig(path + exp + '/' + exp + '_' + scan + '_BFS.png', \
                    transparent=True, bbox_inches='tight', pad_inches=0.01)
        plt.savefig(path + exp + '/PDFs/' + exp + '_' + scan + '_BFS.pdf', \
                    transparent=True, bbox_inches='tight', pad_inches=0.01)

    ### Plot linewidth
    plt.matshow(procData[exp][scan]['LWArr'][:,0,:], cmap=colormap, \
                vmin=lmin, vmax=lmax, origin='lower', interpolation=None, \
                extent=(0, xrange, 0, zrange))
    plt.colorbar(label='Brillouin peak linewidth [GHz]', shrink=0.75)
    plt.grid(b=None)
    ax = plt.gca()
    ax.axis('off')
    ax.autoscale_view
    scalebar = ScaleBar(1, 'um', length_fraction=0.3, location='upper left', \
                        color='w', box_color='None')
    ax.add_artist(scalebar)
    if saveFig:
        plt.savefig(path + exp + '/' + exp + '_' + scan + '_LW.png', \
                    transparent=True, bbox_inches='tight', pad_inches=0.01)
    
    ### Plot signal intensity (integrated Brillouin photons in S/AS peaks)
    plt.matshow(procData[exp][scan]['IntegrPhotonsArr'][:,0,:], cmap=colormap, \
                vmin=pmin, vmax=pmax, origin='lower', interpolation=None, \
                extent=(0, xrange, 0, zrange))
    plt.colorbar(label='Integrated S/AS photons [#]', shrink=0.75)
    plt.grid(b=None)
    ax = plt.gca()
    ax.axis('off')
    ax.autoscale_view
    scalebar = ScaleBar(1, 'um', length_fraction=0.3, location='upper left', \
                        color='w', box_color='None')
    ax.add_artist(scalebar)
    if saveFig:
        plt.savefig(path + exp + '/' + exp + '_' + scan + '_photons.png', \
                    transparent=True, bbox_inches='tight', pad_inches=0.01)
        
    ### Plot Brightfield images:
    # Size of single FLIR camera frame in um
    imageSize = 9000.0/rawData[exp][scan]['attrs']['paramtree/Microscope Camera/Magnification']
    for z in range(procData[exp][scan]['FreqArr'].shape[0]):
        depth = z*measStep[2]

        if 'CMOSImage' in rawData[exp][scan].keys():
            image = np.rot90(rawData[exp][scan]['CMOSImage'][z], -1, (1,0)) # Rotate by -90 deg.
        elif 'BrightfieldImage' in rawData[exp][scan].keys():
            image = np.rot90(rawData[exp][scan]['BrightfieldImage'][z], -1, (1,0)) # Rotate by -90 deg.
                    
        plt.figure()
        plt.grid(b=None)
        plt.title(exp + '/' + scan + r', z = %.1f ${\rm \mu m}$' %depth)
        plt.imshow(image, cmap='Greys_r', vmin=np.amin(image), vmax=np.amax(image), \
                    interpolation='nearest', extent=(0, imageSize, 0, imageSize))
        ax = plt.gca()
        ax.axis('off')
        ax.autoscale_view
        scalebar = ScaleBar(1, 'um', length_fraction=0.2, location='upper left', \
                            color=(1.0, 0.0, 0.0), box_color='None')
        ax.add_artist(scalebar)
        # ax.xaxis.set_ticks_position('bottom')
        if saveFig:
            plt.savefig(path + exp + '/Brightfield/' + exp + '_' + scan + r'_z_%d.png' %z, \
                    transparent=True, bbox_inches='tight', pad_inches=0.01)


#### Plot Brightfield images (only) for single exp and scan
def plotBrightfield(path, rawData, procData, exp, scan, surfZ=0, saveFig=False):
    
    frames = np.array([rawData[exp][scan]['attrs']['paramtree/Frame Number/X'],
                        rawData[exp][scan]['attrs']['paramtree/Frame Number/Y'],
                        rawData[exp][scan]['attrs']['paramtree/Frame Number/Z']])

    # Measured step is based on the coordinates recorded by the Zaber stage
    measStep = np.array([0.0, 0.0, 0.0])
    if frames[0]>1:
        measStep[0] = rawData[exp][scan]['MotorCoords'][1, 0] - rawData[exp][scan]['MotorCoords'][0, 0]
    if frames[1]>1:
        measStep[1] = rawData[exp][scan]['MotorCoords'][frames[0], 1] - rawData[exp][scan]['MotorCoords'][0, 1]
    if frames[2]>1:
        measStep[2] = rawData[exp][scan]['MotorCoords'][frames[0]*frames[1], 2] - rawData[exp][scan]['MotorCoords'][0, 2]

    # Size of single FLIR camera frame in um
    imageSize = 9000.0/rawData[exp][scan]['attrs']['paramtree/Microscope Camera/Magnification']
    
    for z in range(procData[exp][scan]['FreqArr'].shape[0]):
        if surfZ==0:
            depth = z*measStep[2]
        else:
            depth = z*measStep[2] + rawData[exp][scan]['MotorCoords'][0,2] - surfZ
        print('Frame #%d, depth = %.1f um' %(z, depth))
        
        ### Plot Brightfield image:
        if 'CMOSImage' in rawData[exp][scan].keys():
            image = np.rot90(rawData[exp][scan]['CMOSImage'][z], -1, (1,0)) # Rotate by -90 deg.
        elif 'BrightfieldImage' in rawData[exp][scan].keys():
            image = np.rot90(rawData[exp][scan]['BrightfieldImage'][z], -1, (1,0)) # Rotate by -90 deg.
                    
        plt.figure()
        plt.grid(b=None)
        plt.title(exp + '/' + scan + r', z = %.1f ${\rm \mu m}$' %depth)
        plt.imshow(image, cmap='Greys_r', vmin=np.amin(image), vmax=np.amax(image), \
                    interpolation='nearest', extent=(0, imageSize, 0, imageSize))
        ax = plt.gca()
        ax.axis('off')
        ax.autoscale_view
        scalebar = ScaleBar(1, 'um', length_fraction=0.2, location='upper left', \
                            color=(1.0, 0.0, 0.0), box_color='None')
        ax.add_artist(scalebar)
        # ax.xaxis.set_ticks_position('bottom')
        if saveFig:
            plt.savefig(path + exp + '/' + exp + '_' + scan + r'_z_%d.png' %z, \
                    transparent=True, bbox_inches='tight', pad_inches=0.01)


