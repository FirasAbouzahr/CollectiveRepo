# Firas Abouzahr, John Cesar, Kyle Klein

'''
this serves as a header file for analysis related to experimental studies with PET
detectors on flash (highly intense) proton beam therapy.
'''

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit
import scipy as sp
from os import walk

import os
from tqdm import tqdm
plt.style.use('Desktop/ChannelPairStyles.mplstyle')

ln2 = 0.6391

def generateTextY(num,max,subplots=False):
    if max > 50 and subplots==False:
        y = [max - 0.04*(i+1)*max for i in range(0,num+1)]
    elif max <= 50 and subplots == False:
        y = [max - 0.06*(i+1)*max for i in range(0,num+1)]
    elif subplots == True:
        y = [0.96*max - 0.08*(i+1)*max for i in range(0,num+1)]
    return y

def getPaths(directory):
    pathlist = []
    for (dirpath, dirnames, filenames) in walk(directory):
        pathlist.append(filenames)
        break
    return pathlist

# hard-coded for a good reason ; conversion bewteen DAQ channel ID and geometric channel ID does not have any decipherable pattern...
indices = {
      0 : (4,7-7),
      1 : (4,7-6),
      2 : (7,7-5),
      3 : (5,7-7),
      4 : (5,7-4),
      5 : (5,7-5),
      6 : (4,7-4),
      7 : (7,7-7),
      8 : (6,7-6),
      9 : (7,7-4),
      10 : (5,7-6),
      11 : (6,7-4),
      12 : (4,7-5),
      13 : (6,7-5),
      14 : (6,7-7),
      15 : (7,7-6),
      16 : (3,7-7),
      17 : (3,7-6),
      18 : (2,7-7),
      19 : (2,7-6),
      20 : (0,7-7),
      21 : (1,7-7),
      22 : (0,7-6),
      23 : (1,7-6),
      24 : (3,7-5),
      25 : (1,7-5),
      26 : (2,7-5),
      27 : (4,7-3),
      28 : (0,7-5),
      29 : (3,7-4),
      30 : (0,7-4),
      31 : (1,7-4),
      32 : (2,7-4),
      33 : (3,7-3),
      34 : (2,7-3),
      35 : (0,7-3),
      36 : (1,7-3),
      37 : (0,7-2),
      38 : (5,7-3),
      39 : (1,7-2),
      40 : (2,7-2),
      41 : (3,7-2),
      42 : (1,7-1),
      43 : (0,7-1),
      44 : (0,7-0),
      45 : (3,7-1),
      46 : (1,7-0),
      47 : (2,7-1),
      48 : (3,7-0),
      49 : (2,7-0),
      50 : (6,7-2),
      51 : (6,7-1),
      52 : (7,7-1),
      53 : (4,7-1),
      54 : (5,7-1),
      55 : (6,7-0),
      56 : (7,7-0),
      57 : (7,7-2),
      58 : (7,7-3),
      59 : (4,7-2),
      60 : (5,7-0),
      61 : (5,7-2),
      62 : (6,7-3),
      63 : (4,7-0),
      64:(3+8,7),
      65:(3+8,6),
      66:(2+8,4),
      67:(2+8,6),
      68:(3+8,4),
      69:(1+8,7),
      70:(1+8,5),
      71:(0+8,7),
      72:(1+8,6),
      73:(3+8,3),
      74:(2+8,7),
      75:(2+8,3),
      76:(3+8,5),
      77:(0+8,5),
      78:(2+8,5),
      79:(0+8,6),
      80:(4+8,7),
      81:(6+8,7),
      82:(5+8,7),
      83:(7+8,7),
      84:(5+8,6),
      85:(4+8,6),
      86:(6+8,6),
      87:(7+8,6),
      88:(4+8,5),
      89:(6+8,5),
      90:(5+8,5),
      91:(1+8,4),
      92:(7+8,5),
      93:(7+8,4),
      94:(6+8,4),
      95:(4+8,4),
      96:(5+8,4),
      97:(5+8,3),
      98:(6+8,3),
      99:(4+8,3),
      100:(7+8,3),
      101:(7+8,2),
      102:(0+8,4),
      103:(6+8,2),
      104:(7+8,1),
      105:(5+8,2),
      106:(6+8,1),
      107:(4+8,2),
      108:(7+8,0),
      109:(5+8,1),
      110:(6+8,0),
      111:(4+8,1),
      112:(5+8,0),
      113:(4+8,0),
      114:(0+8,2),
      115:(2+8,1),
      116:(0+8,1),
      117:(3+8,1),
      118:(1+8,1),
      119:(1+8,0),
      120:(0+8,0),
      121:(1+8,2),
      122:(1+8,3),
      123:(3+8,2),
      124:(2+8,0),
      125:(2+8,2),
      126:(0+8,3),
      127:(3+8,0)}

def toGeo(x):
    #converts PETSys ID to geometric ID
    y = 8*indices.get(x)[0] + indices.get(x)[1]
    return y

geo_channels = []
for i in range(128):
    geo_channels.append([i,toGeo(i)])
geo_channels = np.asarray(geo_channels)

def toGeoChannelID(AbsChannelID):
    # Convert PETSys absolute channel IDs to geomteric IDs
    slaveID = AbsChannelID // 4096
    chipID = (AbsChannelID - slaveID*4096) // 64
    channelID = AbsChannelID % 64

    PCB_ChanID = 64*(chipID % 2) + channelID
    AbsPCB_ChanID = geo_channels[geo_channels[:,0] == PCB_ChanID][0][1]

    #General formula can be found in above function "to AbsChannelID"
    GeoChannelID = 10**4 * slaveID + 10**2 * chipID + AbsPCB_ChanID % 64
    return GeoChannelID

'''
this displays coincidence hit histograms for an entire proton beam PET run (whole run, beam on only, and after the beam has been turned off). This data is further isolated for analysis in BeamStructure.py

inSpillBinWidth [ms]
postSpillBinWidth [s]
'''
def coincidenceHistograms(dir,datafiles,inSpillBinWidth,postSpillBinWidth,runLengths,savefigs = False):
    
    # this grabs isolated spill data (± 250 ms from spill peak), this is data collected while beam was on. Could be used for SPECT
    PETx_spill_list = []
    PETy_spill_list = []
    
    # this grabs isolated post-spill (after beam shutoff) data. Used to analyze istope activation and for actual PET imaging
    PETx_list = []
    PETy_list = []
    
    index = 0
    for file in datafiles:
        file_name = file
        
        if savefigs == True:
            save_dir = dir+file_name[:-6]+"/"
            try:
                os.mkdir(save_dir)
            except FileExistsError:
                print("The directory for storing figures for this run already exists")

        data = np.genfromtxt(dir+"{}".format(file_name), delimiter="\t", usecols=(2,3,4,7,8,9))
        # Input data column order is TimeL, ChargeL, ChannelIDL, TimeR, ChargeR, ChannelIDR

        data[:,0] = data[:,0] / 1000000000000
        data[:,3] = data[:,3] / 1000000000000

        # full run coincidence hits histogram
        print('The current data is ' + file)
        num_bins = 2 * runLengths[index]
        fig = plt.figure(figsize=(16,9))
        values_time,bins_time,params_time = plt.hist(data[:,0], bins=num_bins)
        plt.xlabel("Time in s")
        plt.ylabel("No. Coincidences per 500 ms")
        plt.yscale('log')
        
        if savefigs == True:
            plt.savefig(save_dir+"Coincidences-over-time")
        plt.show()

        #We can zoom in on the max bin in the above coarse coincidences over time plot and look at the spill
        bin_centers_time = 0.5*(bins_time[1:] + bins_time[:-1])
        max_bin = np.argmax(values_time)
        #getting the time of the max value (presumably the spill)
        time_spill = bin_centers_time[max_bin]
        #only looking plus or minus 250 ms around the spill
        spill_data = data[(data[:,0] >= time_spill-0.25) & (data[:,0] <= time_spill+0.25)]
        spill_and_post = data[(data[:,0] >= time_spill-0.05)]
        #500 ms / 5000 = 100 microseconds
        fig = plt.figure(figsize=(16,9))
        
         
        numBins = int(500 / inSpillBinWidth)
        PETyspill,PETxspill,_= plt.hist(spill_data[:,0],bins=numBins)
        plt.xlabel("Time [s]")
        plt.ylabel("No. Coincidences per " + str(inSpillBinWidth) + " ms")
        plt.title("Coincidences over time: Zoomed around the spill")
        
        if savefigs == True:
            plt.savefig(save_dir+"Coincidences-over-time-spill")
        
        plt.show()
        
        PETx_spill_list.append(PETxspill)
        PETy_spill_list.append(PETyspill)
        
        fig,ax = plt.subplots(figsize=(16,9))
        starttime = time_spill-0.05
        timefor = (runLengths[index] - starttime)
        numBins = int(timefor/postSpillBinWidth)
        
        PETy,PETx,_=plt.hist(spill_and_post[:,0],bins=numBins)
        plt.xlabel("Time [s]")
        plt.ylabel("No. Coincidences per " + str(postSpillBinWidth) + " s")
        plt.title("Coincidences over time: Spill & Post-Spill")
        if savefigs == True:
            plt.savefig(save_dir+"Coincidences-over-time-spill")
        plt.yscale('log')
        plt.show()
        PETx_list.append(PETx)
        PETy_list.append(PETy)
        index += 1
        
        # this grabs
    return PETx_spill_list,PETy_spill_list,PETx_list,PETy_list
    
    
    
# this function plots the energy spectrum for a given data set and fits to the photopeak
# this is for cut data (e.g., in-spill data is removed so prompt gammas are not clouding the data)
def energyResponse(data,label,bins,title='',lims = [0,0,0,0],customguess = False,guess = [1,1,1],display = False):
    fig,ax = plt.subplots(figsize = (10,7))
    binheights,bincenters,_ = ax.hist(data,bins = bins,color = 'C0',alpha = 0.5,label = 'Channel ID ' + str(label))
    centers = np.array([0.5 * (bincenters[i] + bincenters[i+1]) for i in range(len(bincenters) - 1)])
    
    A = np.max(binheights)
    mu = np.where(binheights == A)[0][0]
    mu = bincenters[mu]
    
    if customguess == False:
        guess = [A,mu,1]
    
    p,c = curve_fit(gaussian,centers,binheights,p0 = guess)
    xspace = np.linspace(10,40,1000)
    ax.plot(xspace,gaussian(xspace,*p),c = 'black')
    
    plt.xlabel('Energy in DAQ Units')
    plt.ylabel('Counts')
    plt.title(title)
    plt.legend()
    
    ax.set_xlim(lims[0],lims[1])
    ax.set_ylim(lims[2],lims[3])
    
    xlims = list(ax.get_xlim())
    ylims = list(ax.get_ylim())
    lims = xlims + ylims
    
    
    if display == True:
        plt.show()
    else:
        plt.close()
        
    return p,bins,lims


# this function plots the energy spectrum for a given data set and fits to the photopeak
# this is a more sophisticated photopeak fitter where one must select in the data where to look at for the photopeak. This is required for fitting to energy during the beam spill where prompt gammas pollute our coincidence windows
def noisyEnergyResponse(channelData,label,title='',customlims = False,lims = [0,0,0,0],cuts = [25,35],binNum = 100,sigma_cut = 2.5,display = False): # majorChannel choose 'left' or 'right'
    data = channelData
    
    channelData = channelData[channelData >= cuts[0]]
    channelData = channelData[channelData <= cuts[1]]
    

    bins = np.linspace(min(channelData),max(channelData),binNum)

    counts,binx,_ = plt.hist(channelData,bins = bins)
    plt.close()

    centers = np.array([0.5 * (binx[i] + binx[i+1]) for i in range(len(binx) - 1)])
    
    # get best guesses for gaussian fit parameters
    A = max(counts)
    mu = binx[np.where(counts == A)[0][0]]
    
    fig,ax = plt.subplots(figsize = (10,7))
    
    try:
        p,c = curve_fit(gaussian,centers,counts,p0 = [A,mu,1])
        ppcut = p[1] - sigma_cut*p[2]
        xspace = np.linspace(15,40,1000)
        plt.plot(xspace,gaussian(xspace,*p),c='black')
    
    except:
        print('---- Photopeak Fit Failed ----')
        p = [1,1,1]
        ppcut = 0
    
    binwidth = abs((max(bins)-min(bins))/len(bins))
    binsUpdated = np.arange(-5,45,binwidth)
    plt.hist(data,bins = binsUpdated,label = 'Channel ID ' + str(label),color = 'C0',alpha = 0.5)

    
    if customlims == True:
        plt.xlim(lims[0],lims[1])
        plt.ylim(lims[2],lims[3])
    else:
        maximum = max(data)
        plt.xlim(-5,maximum)
        
    plt.xlabel('Energy in DAQ Units')
    plt.ylabel('Counts')
    plt.title(title)
    plt.legend()
    xlims = list(ax.get_xlim())
    ylims = list(ax.get_ylim())
    lims = xlims + ylims
    if display == True:
        plt.show()
    else:
        plt.close()
        
    return p,binsUpdated,lims

