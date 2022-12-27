'''
Firas Abouzahr

This script is used to do pixel by pixel (e.g., single LOR) performance analysis
as well as material behavior analysis for Geant4 simulated PET scanners\

'''

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.constants import c,h

# for fitting...
def gaussian(x,mu,sig,A):
    return A * np.exp(-((x-mu)/sig)**2)

def exp(x,A,q):
    return A*np.exp(q*x)

# retrives copy number of SiPM pixel
def getCopyNumber(df,CopyNum):
    newFrame = df[df.CopyNum == CopyNum]
    return newFrame

# Gives a histogram of depth of interactions for incident gammas... used to understand attenuation for DOI
def getGammaStopping(gammaframe,opticalframe,copyNum,position,color = 'black'):
    pos = []
    pos_temp = []
    gammaframe = gammaframe[gammaframe.CopyNum == copyNum]
    gammaframe = gammaframe[gammaframe.Object == 'crystal']
    gammaframe = gammaframe.drop_duplicates(["EventID"],keep='last')
    opticalframe = opticalframe[opticalframe.CopyNum == copyNum]
    opticalframe = opticalframe.drop_duplicates(["EventID"],keep='last')
    
    for i,j in zip(gammaframe.Pos,gammaframe.EventID):
        pos_temp.append([i,j])
    
    pos_opt = list(opticalframe.EventID)
    print(pos_opt)
            
    # drop brackets around
    for values in range(len(pos_temp)):
        value = str(pos_temp[values][0])
        value = value[1:]
        
        pos_temp[values][0] = float(value)

    
    for interactions in pos_temp:
        if interactions[1] in pos_opt:
            print(interactions[1])
            pos.append(interactions[0])
        
    pos = abs(np.array(pos)) - position

    fig,ax = plt.subplots()
    y,x,_= ax.hist(pos,bins = 20)
    x = x[:len(y)]
    attenlength = np.where(x <= 10)
    attenlength = attenlength[0][len(attenlength[0])-1] + 1
    plt.title('Gamma Interaction Depth',color=color)
    plt.ylabel('Counts',color=color)
    plt.xlabel("Depth of Interaction (mm)",color=color)
    plt.text(.01, .99,'Mean = ' + str(round(np.mean(pos),2)),va='top', transform=ax.transAxes)
    
    fig,ax = plt.subplots()
    preX = x[:attenlength+1]
    preY = y[:attenlength+1]
    postX = x[attenlength:]
    postY = y[attenlength:]
    ax.step(preX,preY,label = 'Area = ' + str(round(sum(preY),2)))
    ax.step(postX,postY,label = 'Area = ' + str(round(sum(postY),2)))
    plt.title('Gamma Interaction Depth',color=color)
    plt.ylabel('Counts',color=color)
    plt.xlabel("Depth of Interaction (mm)",color=color)
    plt.legend()
    
    return 0
    
# gives the coincidence time resolution between two SiPM pixels (one LOR)
# LORpixels is given as a list of 2 SiPM channels' copy numbers in coincidence
# n chooses the threshold number of scintillation photons considered a detection of a gamma
# photopeakcut is used to choose what events to include in time differences based on a photopeak cut, put 0 if you dont want it
def getCTR(df,LORpixels,bins,bounds,photopeakcut,n=5, color = 'black',customwidth = False, width = np.linspace(0,1,1)):
    fig,ax = plt.subplots()
    # get data organized first
    left = getCopyNumber(df,LORpixels[0])
    right = getCopyNumber(df,LORpixels[1])
    eventIDs = list(set(df.EventID))
    
    # find gamma detections and grab time differences
    timeDiff = []
    for eventNum in eventIDs:
        photoElectric = df[df.EventID == eventNum]
    
        leftDf = left[left.EventID == eventNum]
        leftTots = leftDf.EventID.value_counts().to_numpy()
        tL = np.sort(leftDf.Time)
    
        rightDf = right[right.EventID == eventNum]
        rightTots = rightDf.EventID.value_counts().to_numpy()
        tR = np.sort(rightDf.Time)
    
        if len(tR) >= n and len(tL) >= n and leftTots >= photopeakcut[0] and rightTots >= photopeakcut[1]:
            tL = tL[n-1]
            tR = tR[n-1]
            timeDiff.append(tR - tL)
        if eventNum % 20 == 0:
            print("Running... (Event #: " + str(eventNum) + ")")
            
    timeDiff = np.array(timeDiff)
    
    # plot distribution of time differences and curve fit to gaussian
    if customwidth == False:
        bins = np.linspace(timeDiff.min(),timeDiff.max(),bins)
    else:
        bins = width
        
    y,x,_ = plt.hist(timeDiff,bins = bins, )
    centers = np.array([0.5 * (bins[i] + bins[i+1]) for i in range(len(bins) - 1)])
    guess = [np.mean(timeDiff),np.std(timeDiff),max(y)] # only need if our distributions isn't very Gaussian...
    p,c = curve_fit(gaussian,centers,y,guess)
    xspace = np.linspace(p[0] - 3*p[1],p[0] + 3*p[1],1000) # gives bounds of curve based on fitting parameters
    plt.plot(xspace,gaussian(xspace,*p))
    plt.xlim(bounds[0],bounds[1])
    FWHM = abs(2.355 * p[1])
    
    # add statistics, titles, labels to plot
    plt.title('Timing Differences',color = color)
    plt.xlabel('Time Differences (ns)',color = color)
    plt.ylabel('Counts',color = color)
    plt.text(.01, .99,'Mean = ' + str(round(p[0],2)),va='top', transform=ax.transAxes)
    plt.text(.01, .93,'std = ' + str(round(abs(p[1]),2)),va='top', transform=ax.transAxes)
    plt.text(.01, .87,'CTR = ' + str(round(FWHM,4)),va='top', transform=ax.transAxes)
    plt.text(.01, .81,'Counts = ' + str(round(len(timeDiff),2)),va='top', transform=ax.transAxes)
    plt.xticks(color=color)
    plt.yticks(color=color)
    return bins

# returns energy deposition plot and energy resolution for a given SiPM pixel (based on copy number)
def getEres(df,copyNum,bins,cut,customguess = False,guess = [0,1,1],color = 'black'):
    
    # get data from requested SiPM pixel
    sipmData = getCopyNumber(df,copyNum)
    ene = sipmData.EventID.value_counts().to_numpy()
    
    # plotting and curve fitting
    fig,ax = plt.subplots()
    bins = np.linspace(ene.min(),ene.max(),bins)
    yEne,xEne,_ = ax.hist(ene,bins = bins, )
    centers = np.array([0.5 * (bins[i] + bins[i+1]) for i in range(len(bins) - 1)])
    peaklocation = np.where(yEne == yEne.max())[0][0]
    mu_guess = xEne[peaklocation]
    
    if customguess == False:
        guess = [mu_guess,np.std(ene),yEne.max()]
    
    else:
        guess = guess
    
    p,c = curve_fit(gaussian,centers,yEne,p0 = guess)
    xspace = np.linspace(xEne.min(),xEne.max(),1000)
    ax.plot(xspace,gaussian(xspace,*p),c = 'red')
    FWHM = 2.355 * p[1]
    Eres = FWHM / p[0] * 100
    photopeakcut = p[0] - cut*p[1]
    
    # add titles, statistics, and labels
    plt.title('Energy Spectrum',color = color)
    ax.set_ylabel('Event Counts',color = color)
    ax.set_xlabel('Photon Counts',color = color)
    ax.text(.01, .99,'Mean = ' + str(round(p[0],2)),va='top', transform=ax.transAxes)
    ax.text(.01, .94, 'std = ' + str(round(p[1],2)),va='top', transform=ax.transAxes)
    ax.text(.01, .89,'Eres = ' + str(round(Eres,2)),va='top', transform=ax.transAxes)
    ax.text(.01, .84,'Counts = ' + str(sum(yEne)),va='top', transform=ax.transAxes)
    plt.xticks(color=color)
    plt.yticks(color=color)
    return photopeakcut

# gives scintillation light spectrun of our given scintillant
# mainly used as check for material accuracy
def getEmissionSpectrum(df,bins,material,color = 'black'):
    fig,ax = plt.subplots()
    wl = df.wl.to_numpy()
    y,x,_ = plt.hist(df.wl,bins = bins)
    peakWavelength = np.where(y == max(y))
    peakWavelength = peakWavelength[0][0]
    plt.text(.01, .99,'Peak: ' + str(int(x[peakWavelength])) + ' nm',va='top', transform=ax.transAxes)
    plt.title('Simulated %s Emission Spectrum' % material,color = color)
    plt.xlabel('wavelength (nm)',color = color)
    plt.ylabel('Counts',color = color)
    plt.xticks(color=color)
    plt.yticks(color=color)
    return 0

def getTimeDecay(df,timeconstant,A=1,bins=100,copynum=0,EventID=0,allevents=False,color = 'black'):
    
    fig,ax = plt.subplots()
    
    if allevents == False:
        df = df[df.EventID == EventID]
    
    df = getCopyNumber(df,copynum)

        
    t = df.Time.to_numpy()
    bins = np.linspace(min(df.Time),max(df.Time),bins)
    y,x,_ = ax.hist(t,bins = bins)
    centers = np.array([0.5 * (bins[i] + bins[i+1]) for i in range(len(bins) - 1)])
    guess = [A,timeconstant]
    p,c = curve_fit(exp,x[:len(y)],y,p0 = guess)
    
    xspace = np.linspace(x.min(),x.max(),1000)
    ax.plot(xspace,exp(xspace,p[0],p[1]),linewidth=2,c = 'red',linestyle = 'dashed')
    mean = np.mean(t)
    counts = np.size(t)
    decay = -1/p[1]
    ax.text(.7, .99,'Statistics:',va='top', transform=ax.transAxes, weight='bold')
    ax.text(.7, .94,'Mean = ' + str(round(mean,2)),va='top', transform=ax.transAxes)
    ax.text(.7, .89,'Counts = ' + str(counts),va='top', transform=ax.transAxes)
    ax.text(.7, .84,'$\lambda$ = ' + str(round(decay,3)),va='top', transform=ax.transAxes)
    plt.ylabel("Number of Photons",color=color)
    plt.xlabel("Detection Time (ns)",color=color)
    plt.title('Time Until Detection',color=color)
    plt.xticks(color = color)
    plt.yticks(color = color)
    return 0

def getMeanNumCerenkovPhotons(df,):
    cevents = []
    for i in list(set(df.EventID)):
        tempdf = df[df.EventID == i]
        df_ceren = tempdf[tempdf.Process == 'Cerenkov']
        cevents.append(np.shape(df_ceren)[0])
    
    mean = np.mean(cevents)
    return mean

def getPhotoDetectionEff(df,master):
    totalPhotons = np.shape(master)[0]
    totalDetectedPhotons = np.shape(df)[0]
    PDE = totalDetectedPhotons/totalPhotons
    return PDE


def getTimeDiffs(df,LORpixels,photopeakcut,n=5):
    fig,ax = plt.subplots()
    # get data organized first
    left = getCopyNumber(df,LORpixels[0])
    right = getCopyNumber(df,LORpixels[1])
    eventIDs = list(set(df.EventID))
    
    # find gamma detections and grab time differences
    timeDiff = []
    cerenkovTimeDiff = []
    halfCerenkovTimeDiff = []
    fullCerenkovTimeDiff = []
    
    for eventNum in eventIDs:
        photoElectric = df[df.EventID == eventNum]
    
        leftDf = left[left.EventID == eventNum]
        timeL = leftDf.Time.to_numpy()
        leftTots = leftDf.EventID.value_counts().to_numpy()
        tL = np.sort(leftDf.Time)
        processL = leftDf.Process.to_numpy()
    
        rightDf = right[right.EventID == eventNum]
        rightTots = rightDf.EventID.value_counts().to_numpy()
        timeR = rightDf.Time.to_numpy()
        tR = np.sort(rightDf.Time)
        processR = rightDf.Process.to_numpy()
    
        if len(tR) >= n and len(tL) >= n and leftTots >= photopeakcut[0] and rightTots >= photopeakcut[1]:
            tL = tL[n-1]
            tR = tR[n-1]
            timeDiff.append(tR - tL)
            L = np.where(timeL == tL)
            processL = processL[L[0][0]]
            R = np.where(timeR == tR)
            processR = processR[R[0][0]]
            
            if processR == 'Cerenkov' or processL == 'Cerenkov' and processR != processL:
                halfCerenkovTimeDiff.append(tR - tL)
                cerenkovTimeDiff.append(tR - tL)
            elif processR == 'Cerenkov' and processL == 'Cerenkov':
                fullCerenkovTimeDiff.append(tR - tL)
                cerenkovTimeDiff.append(tR - tL)
            
            
        if eventNum % 20 == 0:
            print("Running... (Event #: " + str(eventNum) + ")")
            
    timeDiff = np.array(timeDiff)
    return timeDiff,cerenkovTimeDiff,halfCerenkovTimeDiff,fullCerenkovTimeDiff
