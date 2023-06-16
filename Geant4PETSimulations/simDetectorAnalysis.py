'''
Firas Abouzahr

This script is used to do pixel by pixel (e.g., single line-of-response) performance analysis
as well as material behavior analysis for Geant4 simulated PET scanners

'''

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.constants import c,h

# for fitting...
def gaussian(x,A,mu,sig):
    return A * np.exp(-((x-mu)/sig)**2)

def exp(x,A,q):
    return A*np.exp(q*x)

# retrives copy number of SiPM pixel
def getCopyNumber(df,CopyNum):
    newFrame = df[df.CopyNum == CopyNum]
    return newFrame

# this is so we can streamline how data is read in
# functions use specific column names so all dataframes must adhere to it so the functions work properly
# read in .txt containing data of all detected photons
def getDetectionDataFrame(file,low_memory = False):
    df = pd.read_csv(file,header = None,low_memory = low_memory)
    df.columns = ['EventID','CopyNum','Time','Process','wavelength']
    return df

# read in .txt containing data of all produced scintillation photons
def getScintilationDataFrane(file,low_memory = False,chunksize = 1):
    df = pd.read_csv(file,header = None,low_memory = low_memory)
    df.columns = ['EventID','Energy','Process']
    return df
    

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
    
# DO NOT USE, THIS IS SHIT PROGRAMMING!
def getCTR(df,LORpixels,bins,bounds,photopeakcut,n=5, color = 'black',customwidth = False, width = np.linspace(0,1,1),title = ''):
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
        if eventNum % 500 == 0:
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
    plt.title(title,color = color)
    plt.xlabel('Time Differences (ns)',color = color)
    plt.ylabel('Counts',color = color)
    plt.text(.01, .99,'Mean = ' + str(round(p[0],2)),va='top', transform=ax.transAxes)
    plt.text(.01, .93,'std = ' + str(round(abs(p[1]),2)),va='top', transform=ax.transAxes)
    plt.text(.01, .87,'CTR = ' + str(round(FWHM,4)),va='top', transform=ax.transAxes)
    plt.text(.01, .81,'Counts = ' + str(round(len(timeDiff),2)),va='top', transform=ax.transAxes)
    plt.xticks(color=color)
    plt.yticks(color=color)
    return bins,FWHM

# plots out energy deposition with normal fit and energy resolution
# returns photopeak cut to cut data compton scattering data from coincidence time resolution
def getEres(df,copyNum,bins,cut,guess = [0,1,1],color = 'black'):
    
    # get data from requested SiPM pixel
    sipmData = getCopyNumber(df,copyNum)
    ene = sipmData.EventID.value_counts().to_numpy()
    
    # plotting and curve fitting
    fig,ax = plt.subplots()
    bins = np.linspace(ene.min(),ene.max(),bins)
    yEne,xEne,_ = ax.hist(ene,bins = bins, )
    centers = np.array([0.5 * (bins[i] + bins[i+1]) for i in range(len(bins) - 1)])

    p,c = curve_fit(gaussian,centers,yEne,p0 = guess)
    xspace = np.linspace(xEne.min(),xEne.max(),1000)
    ax.plot(xspace,gaussian(xspace,*p),c = 'red')
    FWHM = 2.355 * p[2]
    Eres = FWHM / p[1] * 100
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

# plots channel pairs together, do not use to get photopeak cut, only use getEres.
def getPairEres(df,copyNums,bins,guess = [0,1,1]):
    
    # get data from requested SiPM pixel
    pixel0 = getCopyNumber(df,copyNums[0])
    pixel1 = getCopyNumber(df,copyNums[1])
    
    ene0 = pixel0.EventID.value_counts().to_numpy()
    mask0 = ene0 > 10
    ene0 = ene0[mask0]
    
    ene1 = pixel1.EventID.value_counts().to_numpy()
    mask1 = ene1 > 10
    ene1 = ene1[mask1]
    
    # plotting and curve fitting
    fig,ax = plt.subplots()
    bins = np.linspace(ene0.min(),ene0.max(),bins)
    
    yEne0,xEne0,_ = ax.hist(ene0,bins = bins,label = 'Channel ' + str(copyNums[0]) + ' Spectrum',color = 'green',alpha = 0.5)
    yEne1,xEne1,_ = ax.hist(ene1,bins = bins,label = 'Channel ' + str(copyNums[1]) + ' Spectrum',color = 'blue',alpha = 0.5)
    centers = np.array([0.5 * (bins[i] + bins[i+1]) for i in range(len(bins) - 1)])
    
    p0,c0 = curve_fit(gaussian,centers,yEne0,p0 = guess)
    p1,c1 = curve_fit(gaussian,centers,yEne1,p0 = guess)
    
    xspace = np.linspace(xEne0.min(),xEne0.max(),1000)
    ax.plot(xspace,gaussian(xspace,*p0),c = 'red',linestyle = 'dashed')
    ax.plot(xspace,gaussian(xspace,*p1),c = 'red',linestyle = 'dashed')
    FWHM0 = 2.355 * p0[1]
    Eres0 = abs(FWHM0 / p0[0]) * 100
    FWHM1 = 2.355 * p1[1]
    Eres1 = abs(FWHM1 / p1[0]) * 100
    
    # add titles, statistics, and labels
    plt.title('Energy Spectrum (' + str(copyNums[0]) + ',' + str(copyNums[1]) + ')' )
    ax.set_ylabel('Event Counts')
    ax.set_xlabel('Photon Counts')
    
    ax.text(.01, .99,'Stats (Channel '  + str(copyNums[0]) + ')',va='top', transform=ax.transAxes)
    ax.text(.01, .94,'Mean = ' + str(round(p0[0],2)),va='top', transform=ax.transAxes)
    ax.text(.01, .89, 'std = ' + str(round(abs(p0[1]),2)),va='top', transform=ax.transAxes)
    ax.text(.01, .84,'Eres = ' + str(round(Eres0,2)),va='top', transform=ax.transAxes)
    ax.text(.01, .79,'Counts = ' + str(sum(yEne0)),va='top', transform=ax.transAxes)
    ax.text(.01, .71,'Stats (Channel '  + str(copyNums[1]) + ')',va='top', transform=ax.transAxes)
    ax.text(.01, .66,'Mean = ' + str(round(p1[0],2)),va='top', transform=ax.transAxes)
    ax.text(.01, .61, 'std = ' + str(round(abs(p1[1]),2)),va='top', transform=ax.transAxes)
    ax.text(.01, .56,'Eres = ' + str(round(Eres1,2)),va='top', transform=ax.transAxes)
    ax.text(.01, .51,'Counts = ' + str(sum(yEne1)),va='top', transform=ax.transAxes)
    plt.xticks()
    plt.yticks()
    fig.legend(bbox_to_anchor = (1.05,1))
    return 0

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

def getMeanNumCerenkovPhotons(df):
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


# this takes awhile to run but helps eliminate channel pairs connected by compton scatters
def getMostPopulatedChannelPairs(df):
    pairs = []
    eventIDs = list(set(df.EventID))
    for event in eventIDs:
        dfevent = df[df.EventID == event]
        nums = dfevent.CopyNum.to_numpy()
        ch1 = np.bincount(nums).argmax()
        tempnums = np.delete(nums,np.where(nums == ch1))
        if len(tempnums) != 0:
            ch2 = np.bincount(tempnums).argmax()
            if [ch1,ch2] not in pairs and [ch2,ch1] not in pairs:
                pairs.append([ch1,ch2])
        
    return pairs
