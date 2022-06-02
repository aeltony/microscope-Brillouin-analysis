#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 6 15:34:54 2021

@author: ame38
"""
import numpy as np
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
import itertools

#### Fit Brillouin spectrum,
# sline is the data (counts) for the pixels on the spectral line,
# ftol and xtol are fitting tolerances (adjust for speed vs. accuracy)
# rSq_thresh is a minimum r^2 threshold (to assess fit quality)
def fitSpectrum(sline, win_fac=0, rSq_thresh=0.1, xtol=1e-6, ftol=1e-6, verbose=False):
    slineMax = np.amax(np.abs(sline))
    slineCtr = np.floor(0.5*sline.shape[0])
    pix = np.arange(0, sline.shape[0]) # Pixel number
    # Find peak locations
    prominence = 0.4*slineMax  # 0.3*slineMax
    wlen = 5*prominence
    pk_ind, pk_info = find_peaks(sline, prominence=prominence, width=2, \
        height=50, rel_height=0.5, wlen=wlen)
    pk_wids = 0.5*pk_info['widths']
    pk_hts = pk_info['peak_heights']
    
    # Remove non-sensical peaks:
    pk_ind = pk_ind[pk_hts>0]
    pk_wids = pk_wids[pk_hts>0]
    pk_hts = pk_hts[pk_hts>0]
    pk_ind = pk_ind[pk_hts<2*slineMax]
    pk_wids = pk_wids[pk_hts<2*slineMax]
    pk_hts = pk_hts[pk_hts<2*slineMax]
    pk_ind = pk_ind[pk_wids>0]
    pk_wids = pk_wids[pk_wids>0]
    pk_hts = pk_hts[pk_wids>0]
    pk_ind = pk_ind[pk_wids<100]
    pk_wids = pk_wids[pk_wids<100]
    pk_hts = pk_hts[pk_wids<100]

    # Check for no peaks:
    if len(pk_ind)<1:
        if verbose:
            print('[DataFitting] Warning: Too few peaks in spectrum')
        interPeaksteps = np.nan
        linewidth = np.nan
        interPeaksteps_sig = np.nan
        linewidth_sig = np.nan
        fittedSpect = np.nan*np.ones(sline.shape)
        return (interPeaksteps, linewidth, interPeaksteps_sig, linewidth_sig, fittedSpect)

    # Check for single peak (i.e. overlapping peaks):
    elif len(pk_ind)==1:
        # First, check if single peak is centrally located
        pk_ctr_dist = pk_ind - slineCtr
        if np.abs(pk_ctr_dist) < 10:
            # Fit spectrum to 1-Lorentzian model:
            if verbose:
                print('[DataFitting] Warning: Overlapping peaks detected (1-peak fit)')
            # Starting guesses for fit:
            p0 = [pk_hts[0], slineCtr, pk_wids[0], \
                  np.amin(sline)]
            
            # Create boolean mask to filter out points far from the peaks:
            pk_mask = np.array(0*sline, dtype=bool)
            if win_fac > 0:
                pk_mask[(pk_ind[0] - win_fac*pk_wids[0]).astype(int):(pk_ind[0] + win_fac*pk_wids[0]).astype(int)]=True
            else:
                pk_mask[:]=True # Do not use mask

            # Fit spectrum to 1-Lorentzian model:
            try:
                popt, pcov = curve_fit(_1Lorentzian, pix[pk_mask], sline[pk_mask], \
                                  p0=p0, ftol=ftol, xtol=xtol)
                psig = np.sqrt(np.diag(pcov))

                interPeaksteps = 0
                linewidth = np.abs(popt[2])
                interPeaksteps_sig = psig[1]
                linewidth_sig = psig[2]
                fittedSpect = _1Lorentzian(pix, popt[0], popt[1], popt[2], popt[3])
                residuals = sline - fittedSpect
                ss_res = np.sum(residuals**2)
                ss_tot = np.sum((sline-np.mean(sline))**2)
                r_squared = 1 - (ss_res / ss_tot)
            except:
                if verbose:
                    print('[DataFitting] Warning: Fitting 1-peak spectrum failed')
                interPeaksteps = np.nan
                linewidth = np.nan
                interPeaksteps_sig = np.nan
                linewidth_sig = np.nan
                fittedSpect = np.nan*np.ones(sline.shape)
                return (interPeaksteps, linewidth, interPeaksteps_sig, linewidth_sig, fittedSpect)
        
            # R^2 quality check
            if r_squared < rSq_thresh:
                if verbose:
                    print('[DataFitting] Warning: Fitting 1-peak spectrum failed - low R^2')
                interPeaksteps = np.nan
                linewidth = np.nan
                interPeaksteps_sig = np.nan
                linewidth_sig = np.nan
                # fittedSpect = np.nan*np.ones(sline.shape)
            return (interPeaksteps, linewidth, interPeaksteps_sig, linewidth_sig, fittedSpect)

        else:
            if verbose:
                print('[DataFitting] Warning: Potentially overlapping peaks')
            pk_ind = np.array([slineCtr-10, slineCtr+10])  # pk_ind[0]-2
            pk_wids = np.array([pk_wids[0], pk_wids[0]])
            pk_hts = np.array([pk_hts[0], pk_hts[0]])

    # Check for more than 2 peaks:
    elif len(pk_ind)>2:
        if verbose:
            print('[DataFitting] Warning: Potentially > 2 peaks')

        # Remove significantly smaller peaks
        pk_srt = np.argsort(pk_hts)
        pk_select = pk_hts/np.nanmean(pk_hts[pk_srt[-2:]]) > 0.4  # 0.2
        pk_ind = pk_ind[pk_select]
        pk_wids = pk_wids[pk_select]
        pk_hts = pk_hts[pk_select]
        
        # If 3 or 5 peaks remain, remove asymmetric peak
        if len(pk_ind)==3 or len(pk_ind)==5:
            pk_ctr_dist = pk_ind - slineCtr
            pairs = np.array([(a,b) for a,b in itertools.combinations(np.arange(0, len(pk_ind)), 2)])
            sums = np.array([pk_ctr_dist[a]+pk_ctr_dist[b] for (a,b) in pairs])
            sum_srt = np.argsort(np.abs(sums))
            if len(pk_ind)==3:
                # First, check for central overlapping peaks (i.e. 2 pairs of peaks total)
                if np.amin(np.abs(pk_ctr_dist)) < 10:
                    # Keep all 3 peaks for now
                    pk_select = np.unique(pairs)
                else:
                    pk_select = np.unique(pairs[sum_srt[:1]])
            else:
                pk_select = np.unique(pairs[sum_srt[:2]])
                while len(pk_select)<4:
                    sum_srt = np.delete(sum_srt, 1)
                    pk_select = np.unique(pairs[sum_srt[:2]])
            pk_ind = pk_ind[pk_select]
            pk_wids = pk_wids[pk_select]
            pk_hts = pk_hts[pk_select]
    
    # Check for central overlapping peaks + additional pair of peaks
    if len(pk_ind)==3:
        if verbose:
            print('[DataFitting] Warning: Overlapping peaks + 2nd pair of peaks detected')
        # In this case, fit to 3-Lorentzian model:
        # Starting guesses for fit:
        p0 = [pk_hts[0], pk_ind[0], pk_wids[0], \
              pk_hts[1], pk_ind[1], pk_wids[1], \
              pk_hts[2], pk_ind[2], pk_wids[2], \
              np.amin(sline)]
        
        # Create boolean mask to filter out points far from the peaks:
        pk_mask = np.array(0*sline, dtype=bool)
        if win_fac > 0:
            pk_mask[(pk_ind[0] - win_fac*pk_wids[0]).astype(int):(pk_ind[0] + win_fac*pk_wids[0]).astype(int)]=True
            pk_mask[(pk_ind[1] - win_fac*pk_wids[1]).astype(int):(pk_ind[1] + win_fac*pk_wids[1]).astype(int)]=True
            pk_mask[(pk_ind[2] - win_fac*pk_wids[2]).astype(int):(pk_ind[2] + win_fac*pk_wids[2]).astype(int)]=True
        else:
            pk_mask[:]=True # Do not use mask
        
        # Fit spectrum to 3-Lorentzian model
        try:
            popt, pcov = curve_fit(_3Lorentzian, pix[pk_mask], sline[pk_mask], \
                              p0=p0, ftol=ftol, xtol=xtol)
            psig = np.sqrt(np.diag(pcov))
            
            # Remove any non-sensical peaks
            pk_ht_ind = np.array([0, 3, 6])
            pk_ht_ind = pk_ht_ind[popt[pk_ht_ind] > 0]
            pk_ht_ind = pk_ht_ind[popt[pk_ht_ind] < 2*slineMax]
            pk_ht_ind = pk_ht_ind[np.abs(popt[pk_ht_ind+2]) < 100]
            
            # Check if 2 peaks remain
            if pk_ht_ind.shape[0] == 2:
                 pk_1_ind = pk_ht_ind[0]
                 pk_2_ind = pk_ht_ind[1]
                 interPeaksteps = np.abs(popt[pk_1_ind+1] - popt[pk_2_ind+1])
                 linewidth = 0.5*(np.abs(popt[pk_1_ind+2]) + np.abs(popt[pk_2_ind+2])) # Mean linewidth of 2 peaks
                 interPeaksteps_sig = np.sqrt(psig[pk_1_ind+1]**2 + psig[pk_2_ind+1]**2)
                 linewidth_sig = 0.5*np.sqrt(psig[pk_1_ind+2]**2 + psig[pk_2_ind+2]**2)
                 fittedSpect = _2Lorentzian(pix, popt[pk_1_ind], popt[pk_1_ind+1], popt[pk_1_ind+2], \
                                            popt[pk_2_ind], popt[pk_2_ind+1], popt[pk_2_ind+2], popt[9])
                 residuals = sline - fittedSpect
                 ss_res = np.sum(residuals**2)
                 ss_tot = np.sum((sline-np.mean(sline))**2)
                 r_squared = 1 - (ss_res / ss_tot)  
            # Check if a single peak remains
            elif pk_ht_ind.shape[0] == 1:
                if verbose:
                    print('[DataFitting] Warning: Overlapping peaks detected (3-peak fit - A)')
                interPeaksteps = 0
                linewidth = np.abs(popt[pk_ht_ind+2])
                interPeaksteps_sig = psig[pk_ht_ind+1]
                linewidth_sig = psig[pk_ht_ind+2]
                fittedSpect = _1Lorentzian(pix, popt[pk_ht_ind], popt[pk_ht_ind+1], popt[pk_ht_ind+2], popt[9])
                residuals = sline - fittedSpect
                ss_res = np.sum(residuals**2)
                ss_tot = np.sum((sline-np.mean(sline))**2)
                r_squared = 1 - (ss_res / ss_tot)
            else:
                # Find larger pair of peaks
                pk_ctr_dist = np.abs(popt[pk_ht_ind+1] - slineCtr)
                pk_dist_srt = np.argsort(pk_ctr_dist)
                if pk_ht_ind[pk_dist_srt[0]] == pk_ht_ind[np.argmax(popt[pk_ht_ind])]:
                    # Select only the more central (overlapped) peak
                    if verbose:
                        print('[DataFitting] Warning: Overlapping peaks detected (3-peak fit - B)')
                    overlapped_pk_ind = pk_ht_ind[pk_dist_srt[0]]
                    interPeaksteps = 0
                    linewidth = np.abs(popt[overlapped_pk_ind+2])
                    interPeaksteps_sig = psig[overlapped_pk_ind+1]
                    linewidth_sig = psig[overlapped_pk_ind+2]
                    fittedSpect = _1Lorentzian(pix, popt[overlapped_pk_ind], popt[overlapped_pk_ind+1], popt[overlapped_pk_ind+2], \
                                                popt[9])
                    # fittedSpect = _3Lorentzian(pix[pk_mask], \
                    #                 popt[0], popt[1], popt[2], popt[3], popt[4], \
                    #                 popt[5], popt[6], popt[7], popt[8], popt[9])
                    residuals = sline - fittedSpect
                    ss_res = np.sum(residuals**2)
                    ss_tot = np.sum((sline-np.mean(sline))**2)
                    r_squared = 1 - (ss_res / ss_tot)
                else:
                    # Select outer pair of (non-overlapped) peaks
                    pk_1_ind = pk_ht_ind[pk_dist_srt[-1]]
                    pk_2_ind = pk_ht_ind[pk_dist_srt[-2]]
                    interPeaksteps = np.abs(popt[pk_1_ind+1] - popt[pk_2_ind+1])
                    linewidth = 0.5*(np.abs(popt[pk_1_ind+2]) + np.abs(popt[pk_2_ind+2])) # Mean linewidth of 2 peaks
                    interPeaksteps_sig = np.sqrt(psig[pk_1_ind+1]**2 + psig[pk_2_ind+1]**2)
                    linewidth_sig = 0.5*np.sqrt(psig[pk_1_ind+2]**2 + psig[pk_2_ind+2]**2)
                    fittedSpect = _2Lorentzian(pix, popt[pk_1_ind], popt[pk_1_ind+1], popt[pk_1_ind+2], \
                                                popt[pk_2_ind], popt[pk_2_ind+1], popt[pk_2_ind+2], popt[9])
                    # fittedSpect = _3Lorentzian(pix[pk_mask], \
                    #                 popt[0], popt[1], popt[2], popt[3], popt[4], \
                    #                 popt[5], popt[6], popt[7], popt[8], popt[9])
                    residuals = sline - fittedSpect
                    ss_res = np.sum(residuals**2)
                    ss_tot = np.sum((sline-np.mean(sline))**2)
                    r_squared = 1 - (ss_res / ss_tot)

        except:
            if verbose:
                print('[DataFitting] Warning: Fitting 3-peak spectrum failed')
            interPeaksteps = np.nan
            linewidth = np.nan
            interPeaksteps_sig = np.nan
            linewidth_sig = np.nan
            fittedSpect = np.nan*np.ones(sline.shape)
            return (interPeaksteps, linewidth, interPeaksteps_sig, linewidth_sig, fittedSpect)
        
        # R^2 quality check
        if r_squared < rSq_thresh:
            if verbose:
                print('[DataFitting] Warning: Fitting 3-peak spectrum failed - low R^2')
            interPeaksteps = np.nan
            linewidth = np.nan
            interPeaksteps_sig = np.nan
            linewidth_sig = np.nan
            # fittedSpect = np.nan*np.ones(sline.shape)
        return (interPeaksteps, linewidth, interPeaksteps_sig, linewidth_sig, fittedSpect)

    # Check for 2 pairs of peaks
    elif len(pk_ind)==4:
        if verbose:
            print('[DataFitting] Warning: Two pairs of peaks detected')
        ### Fit spectrum to 4-Lorentzian model:
        # Starting guesses for fit:
        p0 = [pk_hts[0], pk_ind[0], pk_wids[0], \
              pk_hts[1], pk_ind[1], pk_wids[1], \
              pk_hts[2], pk_ind[2], pk_wids[2], \
              pk_hts[3], pk_ind[3], pk_wids[3], \
              np.amin(sline)]
        
        # Create boolean mask to filter out points far from the peaks:
        pk_mask = np.array(0*sline, dtype=bool)
        if win_fac > 0:
            pk_mask[(pk_ind[0] - win_fac*pk_wids[0]).astype(int):(pk_ind[0] + win_fac*pk_wids[0]).astype(int)]=True
            pk_mask[(pk_ind[1] - win_fac*pk_wids[1]).astype(int):(pk_ind[1] + win_fac*pk_wids[1]).astype(int)]=True
            pk_mask[(pk_ind[2] - win_fac*pk_wids[2]).astype(int):(pk_ind[2] + win_fac*pk_wids[2]).astype(int)]=True
            pk_mask[(pk_ind[3] - win_fac*pk_wids[3]).astype(int):(pk_ind[3] + win_fac*pk_wids[3]).astype(int)]=True
        else:
            pk_mask[:]=True # Do not use mask
        
        # Fit spectrum to 4-Lorentzian model
        try:
            popt, pcov = curve_fit(_4Lorentzian, pix[pk_mask], sline[pk_mask], \
                              p0=p0, ftol=ftol, xtol=xtol)
            psig = np.sqrt(np.diag(pcov))
            
            # Remove any non-sensical peaks
            pk_ht_ind = np.array([0, 3, 6, 9])
            pk_ht_ind = pk_ht_ind[popt[pk_ht_ind] > 0]
            pk_ht_ind = pk_ht_ind[popt[pk_ht_ind] < 2*slineMax]
            pk_ht_ind = pk_ht_ind[np.abs(popt[pk_ht_ind+2]) < 100]
            
            # Find larger pair of peaks
            pk_ctr_dist = np.abs(popt[pk_ht_ind+1] - slineCtr)
            pk_dist_srt = np.argsort(pk_ctr_dist)
            if np.mean(popt[pk_ht_ind[pk_dist_srt[-2:]]]) > np.mean(popt[pk_ht_ind[pk_dist_srt[:2]]]):
                pk_1_ind = pk_ht_ind[pk_dist_srt[-1]]
                pk_2_ind = pk_ht_ind[pk_dist_srt[-2]]
                interPeaksteps = np.abs(popt[pk_1_ind+1] - popt[pk_2_ind+1])
                linewidth = 0.5*(np.abs(popt[pk_1_ind+2]) + np.abs(popt[pk_2_ind+2])) # Mean linewidth of 2 peaks
                interPeaksteps_sig = np.sqrt(psig[pk_1_ind+1]**2 + psig[pk_2_ind+1]**2)
                linewidth_sig = 0.5*np.sqrt(psig[pk_1_ind+2]**2 + psig[pk_2_ind+2]**2)
                fittedSpect = _2Lorentzian(pix, popt[pk_1_ind], popt[pk_1_ind+1], popt[pk_1_ind+2], \
                                            popt[pk_2_ind], popt[pk_2_ind+1], popt[pk_2_ind+2], popt[12])
                # fittedSpect = _4Lorentzian(pix[pk_mask], \
                #                 popt[0], popt[1], popt[2], popt[3], popt[4], \
                #                 popt[5], popt[6], popt[7], popt[8], popt[9], \
                #                 popt[10], popt[11], popt[12])
                residuals = sline - fittedSpect
                ss_res = np.sum(residuals**2)
                ss_tot = np.sum((sline-np.mean(sline))**2)
                r_squared = 1 - (ss_res / ss_tot)
            else:
                pk_1_ind = pk_ht_ind[pk_dist_srt[0]]
                pk_2_ind = pk_ht_ind[pk_dist_srt[1]]
                interPeaksteps = np.abs(popt[pk_1_ind+1] - popt[pk_2_ind+1])
                linewidth = 0.5*(np.abs(popt[pk_1_ind+2]) + np.abs(popt[pk_2_ind+2])) # Mean linewidth of 2 peaks
                interPeaksteps_sig = np.sqrt(psig[pk_1_ind+1]**2 + psig[pk_2_ind+1]**2)
                linewidth_sig = 0.5*np.sqrt(psig[pk_1_ind+2]**2 + psig[pk_2_ind+2]**2)
                fittedSpect = _2Lorentzian(pix, popt[pk_1_ind], popt[pk_1_ind+1], popt[pk_1_ind+2], \
                                            popt[pk_2_ind], popt[pk_2_ind+1], popt[pk_2_ind+2], popt[12])
                # fittedSpect = _4Lorentzian(pix[pk_mask], \
                #                 popt[0], popt[1], popt[2], popt[3], popt[4], \
                #                 popt[5], popt[6], popt[7], popt[8], popt[9], \
                #                 popt[10], popt[11], popt[12])
                residuals = sline - fittedSpect
                ss_res = np.sum(residuals**2)
                ss_tot = np.sum((sline-np.mean(sline))**2)
                r_squared = 1 - (ss_res / ss_tot)
            
            # # Find larger pair of peaks
            # pk_ht_srt = np.argsort(popt[pk_ht_ind])
            # pk_1_ind = pk_ht_ind[pk_ht_srt[-1]]
            # pk_2_ind = pk_ht_ind[pk_ht_srt[-2]]
            # interPeaksteps = np.abs(popt[pk_1_ind+1] - popt[pk_2_ind+1])
            # linewidth = 0.5*(np.abs(popt[pk_1_ind+2]) + np.abs(popt[pk_2_ind+2])) # Mean linewidth of 2 peaks
            # interPeaksteps_sig = np.sqrt(psig[pk_1_ind+1]**2 + psig[pk_2_ind+1]**2)
            # linewidth_sig = 0.5*np.sqrt(psig[pk_1_ind+2]**2 + psig[pk_2_ind+2]**2)
            # # fittedSpect = _2Lorentzian(pix, popt[pk_1_ind], popt[pk_1_ind+1], popt[pk_1_ind+2], \
            # #                             popt[pk_2_ind], popt[pk_2_ind+1], popt[pk_2_ind+2], popt[12])
            # fittedSpect = _4Lorentzian(pix[pk_mask], \
            #                 popt[0], popt[1], popt[2], popt[3], popt[4], \
            #                 popt[5], popt[6], popt[7], popt[8], popt[9], \
            #                 popt[10], popt[11], popt[12])
            # residuals = sline - fittedSpect
            # ss_res = np.sum(residuals**2)
            # ss_tot = np.sum((sline-np.mean(sline))**2)
            # r_squared = 1 - (ss_res / ss_tot)
            
        except:
            if verbose:
                print('[DataFitting] Warning: Fitting 4-peak spectrum failed')
            interPeaksteps = np.nan
            linewidth = np.nan
            interPeaksteps_sig = np.nan
            linewidth_sig = np.nan
            fittedSpect = np.nan*np.ones(sline.shape)
            return (interPeaksteps, linewidth, interPeaksteps_sig, linewidth_sig, fittedSpect)
        
        # R^2 quality check
        if r_squared < rSq_thresh:
            if verbose:
                print('[DataFitting] Warning: Fitting 4-peak spectrum failed - low R^2')
            interPeaksteps = np.nan
            linewidth = np.nan
            interPeaksteps_sig = np.nan
            linewidth_sig = np.nan
            # fittedSpect = np.nan*np.ones(sline.shape)
        return (interPeaksteps, linewidth, interPeaksteps_sig, linewidth_sig, fittedSpect)
                
    else:
        ### Fit spectrum to 2-Lorentzian model:
        # Starting guesses for fit:
        p0 = [pk_hts[0], pk_ind[0], pk_wids[0], \
              pk_hts[1], pk_ind[1], pk_wids[1], \
              np.amin(sline)]
        
        # Create boolean mask to filter out points far from the peaks:
        pk_mask = np.array(0*sline, dtype=bool)
        if win_fac > 0:
            pk_mask[(pk_ind[0] - win_fac*pk_wids[0]).astype(int):(pk_ind[0] + win_fac*pk_wids[0]).astype(int)]=True
            pk_mask[(pk_ind[1] - win_fac*pk_wids[1]).astype(int):(pk_ind[1] + win_fac*pk_wids[1]).astype(int)]=True
        else:
            pk_mask[:]=True # Do not use mask
    
        # Fit spectrum to 2-Lorentzian model:
        try:
            popt, pcov = curve_fit(_2Lorentzian, pix[pk_mask], sline[pk_mask], \
                              p0=p0, ftol=ftol, xtol=xtol)
            psig = np.sqrt(np.diag(pcov))
            
            # Remove any non-sensical peaks
            pk_ht_ind = np.array([0, 3])
            pk_ht_ind = pk_ht_ind[popt[pk_ht_ind] > 0]
            pk_ht_ind = pk_ht_ind[popt[pk_ht_ind] < 2*slineMax]
            pk_ht_ind = pk_ht_ind[np.abs(popt[pk_ht_ind+2]) < 100]
            
            pk_ctr_dist = popt[pk_ht_ind + 1] - slineCtr
            
            # Check if a single peak remains
            if pk_ht_ind.shape[0] == 1:
                if verbose:
                    print('[DataFitting] Warning: Overlapping peaks detected (2-peak fit - A)')
                interPeaksteps = 0
                linewidth = np.abs(popt[pk_ht_ind+2])
                interPeaksteps_sig = psig[pk_ht_ind+1]
                linewidth_sig = psig[pk_ht_ind+2]
                fittedSpect = _1Lorentzian(pix, popt[pk_ht_ind], popt[pk_ht_ind+1], popt[pk_ht_ind+2], popt[6])
                residuals = sline - fittedSpect
                ss_res = np.sum(residuals**2)
                ss_tot = np.sum((sline-np.mean(sline))**2)
                r_squared = 1 - (ss_res / ss_tot)
            
            # Check for asymmetric peaks (indicating overlapping peaks + a noise peak)
            elif np.abs(pk_ctr_dist[0] + pk_ctr_dist[1]) > 15:
                # Select only the more central (overlapped) peak
                if verbose:
                    print('[DataFitting] Warning: Overlapping peaks detected (2-peak fit - B)')
                if np.abs(pk_ctr_dist[1]) > np.abs(pk_ctr_dist[0]):
                    interPeaksteps = 0
                    linewidth = np.abs(popt[2])
                    interPeaksteps_sig = psig[1]
                    linewidth_sig = psig[2]
                    fittedSpect = _1Lorentzian(pix, popt[0], popt[1], popt[2], popt[6])
                    residuals = sline - fittedSpect
                    ss_res = np.sum(residuals**2)
                    ss_tot = np.sum((sline-np.mean(sline))**2)
                    r_squared = 1 - (ss_res / ss_tot)
                else:
                    interPeaksteps = 0
                    linewidth = np.abs(popt[5])
                    interPeaksteps_sig = psig[4]
                    linewidth_sig = psig[5]
                    fittedSpect = _1Lorentzian(pix, popt[3], popt[4], popt[5], popt[6])
                    residuals = sline - fittedSpect
                    ss_res = np.sum(residuals**2)
                    ss_tot = np.sum((sline-np.mean(sline))**2)
                    r_squared = 1 - (ss_res / ss_tot)
    
            else:
                interPeaksteps = np.abs(popt[4] - popt[1])
                linewidth = 0.5*(np.abs(popt[2]) + np.abs(popt[5])) # Mean linewidth of 2 peaks
                interPeaksteps_sig = np.sqrt(psig[4]**2 + psig[1]**2)
                linewidth_sig = 0.5*np.sqrt(psig[2]**2 + psig[5]**2)
                fittedSpect = _2Lorentzian(pix, popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6])
                residuals = sline - fittedSpect
                ss_res = np.sum(residuals**2)
                ss_tot = np.sum((sline-np.mean(sline))**2)
                r_squared = 1 - (ss_res / ss_tot)
                
        except:
            if verbose:
                print('[DataFitting] Warning: Fitting 2-peak spectrum failed')
            interPeaksteps = np.nan
            linewidth = np.nan
            interPeaksteps_sig = np.nan
            linewidth_sig = np.nan
            fittedSpect = np.nan*np.ones(sline.shape)
            return (interPeaksteps, linewidth, interPeaksteps_sig, linewidth_sig, fittedSpect)
        
        # R^2 quality check
        if r_squared < rSq_thresh:
            if verbose:
                print('[DataFitting] Warning: Fitting 2-peak spectrum failed - low R^2')
            interPeaksteps = np.nan
            linewidth = np.nan
            interPeaksteps_sig = np.nan
            linewidth_sig = np.nan
            # fittedSpect = np.nan*np.ones(sline.shape)
        return (interPeaksteps, linewidth, interPeaksteps_sig, linewidth_sig, fittedSpect)


#### Fit calibration curve to determine SD and FSR
def fitCalCurve(pxDist, freq, xtol=1e-6, ftol=1e-6, verbose=False):
    # Starting guesses for fit:
    p0 = [0.127, 21.5]
    # Remove NaNs before fitting
    freq = freq[~np.isnan(pxDist)]
    pxDist = pxDist[~np.isnan(pxDist)]
    try:
        popt, pcov = curve_fit(_Linear, pxDist, freq, p0=p0, ftol=ftol, xtol=xtol)
    except:
        if verbose:
            print('[DataFitting] Warning: Fitting calibration curve failed')
        SD = np.nan
        FSR = np.nan
        SD_sig = np.nan
        FSR_sig = np.nan
        return (SD, FSR, SD_sig, FSR_sig)
    residuals = freq - _Linear(pxDist, popt[0], popt[1])
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((freq-np.mean(freq))**2)
    r_squared = 1 - (ss_res / ss_tot)
    # print('[DataFitting] Calibration curve fitting r^2 = %.4f' %r_squared)
    if r_squared < 0.9:
        if verbose:
            print('[DataFitting] Warning: Fitting calibration curve failed - low R^2')
        SD = np.nan
        FSR = np.nan
        SD_sig = np.nan
        FSR_sig = np.nan
        return (SD, FSR, SD_sig, FSR_sig)
    else:
        psig = np.sqrt(np.diag(pcov))
        SD = popt[0]
        FSR = popt[1]
        SD_sig = psig[0]
        FSR_sig = psig[1]
    return (SD, FSR, SD_sig, FSR_sig)


#### Fitting functions
def _1Lorentzian(x, amp1,cen1,wid1, offs):
    return (amp1*wid1**2/((x-cen1)**2+wid1**2)) \
            + offs

def _2Lorentzian(x, amp1,cen1,wid1, amp2,cen2,wid2, offs):
    return (amp1*wid1**2/((x-cen1)**2+wid1**2)) \
            + (amp2*wid2**2/((x-cen2)**2+wid2**2)) \
            + offs

def _3Lorentzian(x, amp1,cen1,wid1, amp2,cen2,wid2, amp3,cen3,wid3, offs):
    return (amp1*wid1**2/((x-cen1)**2+wid1**2)) \
            + (amp2*wid2**2/((x-cen2)**2+wid2**2)) \
            + (amp3*wid3**2/((x-cen3)**2+wid3**2)) \
            + offs

def _4Lorentzian(x, amp1,cen1,wid1, amp2,cen2,wid2, amp3,cen3,wid3, amp4,cen4,wid4, offs):
    return (amp1*wid1**2/((x-cen1)**2+wid1**2)) \
            + (amp2*wid2**2/((x-cen2)**2+wid2**2)) \
            + (amp3*wid3**2/((x-cen3)**2+wid3**2)) \
            + (amp4*wid4**2/((x-cen4)**2+wid4**2)) \
            + offs

def _Linear(x, sd, fsr):
    return 0.5*fsr - 0.5*sd*x


