#Kyle Klein, John Cesar, Firas Abouzahr
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit
import scipy as sp
from iminuit import cost, Minuit
import iminuit as im
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

def toGeo(x):
    #converts PETSys ID to geometric ID
    y = 8*indices.get(x)[0] + indices.get(x)[1]
    return y

def getCongregateChannelData(data):
    q = data[0]
    for i in range(1,len(data)):
        q = np.concatenate((q,data[i]))
    
    return q

def gaussian(x,A,mu,sig):
    return A * np.exp(-((x-mu)/sig)**2)
    
# fit energy spectrum
# cuts is part of my fitting algoritim (isolates photopeak for fitting). The predifined cuts are characteristic of our PET detectors response
# sigma_cut specifies the # of std away from the photopeak you want to cut, this photopeak cut is then return --> important for timing resolution, dosimetry etc
def fitAggregateSpectra(channelData,label,cuts = [25,35],binNum = 100,sigma_cut = 2.5,display = False):
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
    
    p,c = curve_fit(gaussian,centers,counts,p0 = [A,mu,1])
    ppcut = p[1] - sigma_cut*p[2]
    xspace = np.linspace(0,max(data),1000)
    
    fig,ax = plt.subplots(figsize = (16,9))
    binwidth = abs((max(bins)-min(bins))/len(bins))
    
    binsUpdated = np.arange(-5,max(data),binwidth)
        
    plt.hist(data,bins = binsUpdated,label = 'Channel ID ' + str(label),color = 'C0',alpha = 0.5)
    plt.plot(xspace,gaussian(xspace,*p),c='black')
    
    maximum = max(data)
    plt.xlim(-5,maximum)
    plt.xlabel('Energy in DAQ Units')
    plt.ylabel('Counts')
    plt.legend()
    if display == True:
        plt.show()
    else:
        plt.close()
        
    return p,binsUpdated,ppcut
    

# p & bins will be retured from fitAggregateSpectra
def getEnergyCounts(channelData,p,bins):
    y,x,_ = plt.hist(channelData,bins = bins)
    plt.close()
     
    totalCounts = np.sum(y)
    
    PhotoPeakLimits = [p[1]-2.5*p[2],p[1]+2.5*p[2]]
    
    PhotoPeakLeftBound = np.where(x >= PhotoPeakLimits[0])[0][0]
    PhotoPeakRightBound = np.where(x >= PhotoPeakLimits[1])[0][0]
    
    photoPeak = y[PhotoPeakLeftBound:PhotoPeakRightBound+1]
    photoPeakCounts = np.sum(photoPeak)
    
    promptGammas = y[PhotoPeakRightBound:]
    promptGammaCounts = np.sum(promptGammas)
    
    return totalCounts,photoPeakCounts,promptGammaCounts


def getPPcuts(data,p):
    return 

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


##############################################################################################
##############################################################################################
#################################### Analysis code ###########################################
##############################################################################################
##############################################################################################

filenames = ['MDA-PMMA1-3min-210mm-8cmCu_coinc.dat',
             'MDA-PMMA1-3min-210mm-4cmCu_coinc.dat',
             'MDA-PMMA1-3min-210mm-NoCu_coinc.dat',
             'MDA-PMMA2-20min-210mm-NoCu_coinc.dat',
             'MDA-HDPE1-20min-210mm-NoCu_coinc.dat',
             'MDA-PMMA3-15min-150mm-NoCu_coinc.dat',
             'MDA-PMMA4-15min-110mm-NoCu_coinc.dat',
             'MDA-H2O1-15min-210mm-NoCu_coinc.dat',
             'MDA-NoPhantom-5min-210mm-NoCu_coinc.dat',
             'MDA-PMMA5-20min-210mm-NoCu_coinc.dat',
             'MDA-Ni-10min-210mm-NoCu_coinc.dat',
             'MDA-Cu-10min-210mm-NoCu_coinc.dat',
             'MDA-PMMA3-15min-150mm-NoCu-2_coinc.dat']

# this data is grabbed using ProtonBeamSpillStructure.py
spillTimes = [[0,.1017],[7.65663546955644, 7.758135469556439],
 [21.6824066438293, 21.7840066438293],
 [62.608382561041, 62.710382561040994],
 [10.5593494114757, 10.661349411475701],
 [68.50224351489645, 68.60394351489644],
 [52.371305190264536, 52.472905190264534],
 [13.31026374859568, 13.412063748595681],
 [69.17326815942364, 69.27506815942364],
 [23.0679552989942, 23.1695552989942],
 [39.92413222839843, 40.025732228398425],
 [76.49712322660356, 76.59892322660356],
 [100.57056977295451, 100.67186977295451]]

runLengths = [180,180,180,1200,1200,900,900,900,300,1200,600,600,900]

# use this to look at one run at a time
runNum = 7
filenames = [filenames[runNum - 1]] 

dir = "Flash_Therapy/PET_3-5-23/"
index = 0
for file in filenames:
    file_name = file
    
    #For showing ricardo the data run 1 was MidV1, run 2 was MidV2, and run3 was MaxV1
    data = np.genfromtxt(dir+"{}".format(file_name), delimiter="\t", usecols=(2,3,4,7,8,9))
    # Input data column order is TimeL, ChargeL, ChannelIDL, TimeR, ChargeR, ChannelIDR

    data[:,0] = data[:,0] / 1000000000000
    data[:,3] = data[:,3] / 1000000000000
    
    spillStart = spillTimes[runNum - 1][0]
    spillStart = np.where(data[:,0] >= spillStart)[0][0]
    
    spillEnd = spillTimes[runNum - 1][1]
    spillEnd = np.where(data[:,0] >= spillEnd)[0][0]

    #First we must determine when the first spill ocurred
    print('The current data is ' + file,'\n')

    # Count the number of unique channel pairs and keep only those with more than 300  events
    CPIDs = data[:,[2,5]] # Take a slice of "data" containing just channel pairs (ChannelIDL and ChannelIDR)
    UniqueCPs = np.unique(CPIDs, axis=0, return_counts=True) # Creates new array of all the unique channel pairs
    CP_num = len(UniqueCPs[1]) # Number of unique channel pair
    UniqueCPs = np.hstack((UniqueCPs[0],UniqueCPs[1].reshape(CP_num, 1)))
    #Sort the unique CPs by the number of coincidences in them
    UniqueCPs = UniqueCPs[UniqueCPs[:,2].argsort()]
    #reverse the ordering so the most populated CPs are at the begining of the list. God I love numpy indexing
    UniqueCPs = UniqueCPs[::-1]

    mostpopulated_L = UniqueCPs[0][0] 
    mostpopulated_R = UniqueCPs[0][1]
    
    totalL = []
    totalR = []
    totalL_inSpill = []
    totalR_inSpill = []
    totalL_postSpill = []
    totalR_postSpill = []
    CP_times = []

    data_spill = data[spillStart:spillEnd+1]
    data_post = data[spillEnd+1:]
    
    print('Getting the aggregate spectra for the most active left channel pair:')
    for i in tqdm(range(len(UniqueCPs[:,1]))):
        # Define channel_pair array as subset of data matching the UniqueCPs array
            LChannel = UniqueCPs[i][0]
            if LChannel == mostpopulated_L:
                RChannel = UniqueCPs[i][1]
                
                channel_pair = data[data[:,2] == UniqueCPs[i][0]]
                channel_pair_spill = data_spill[data_spill[:,2] == UniqueCPs[i][0]]
                channel_pair_post = data_post[data_post[:,2] == UniqueCPs[i][0]]
                
                channel_pair = channel_pair[channel_pair[:,5] == UniqueCPs[i][1]]
                channel_pair_spill = channel_pair_spill[channel_pair_spill[:,5] == UniqueCPs[i][1]]
                channel_pair_post = channel_pair_post[channel_pair_post[:,5] == UniqueCPs[i][1]]
                
                CP_times.append(channel_pair[channel_pair[:,0].argsort()][0][0])
#                 PlotEnergySpectrum(channel_pair, LChannel, RChannel,display=False)
                totalL.append(channel_pair)
                totalL_inSpill.append(channel_pair_spill)
                totalL_postSpill.append(channel_pair_post)
            
    
    print('Getting the aggregate spectra for the most active right channel pair:')          
    for i in tqdm(range(len(UniqueCPs[:,1]))):
        # Define channel_pair array as subset of data matching the UniqueCPs array
            RChannel = UniqueCPs[i][1]
            if RChannel == mostpopulated_R:
                LChannel = UniqueCPs[i][0]
                
                channel_pair = data[data[:,2] == UniqueCPs[i][0]]
                channel_pair_spill = data_spill[data_spill[:,2] == UniqueCPs[i][0]]
                channel_pair_post = data_post[data_post[:,2] == UniqueCPs[i][0]]
                
                channel_pair = channel_pair[channel_pair[:,5] == UniqueCPs[i][1]]
                channel_pair_spill = channel_pair_spill[channel_pair_spill[:,5] == UniqueCPs[i][1]]
                channel_pair_post = channel_pair_post[channel_pair_post[:,5] == UniqueCPs[i][1]]
                
                CP_times.append(channel_pair[channel_pair[:,0].argsort()][0][0])
                totalR.append(channel_pair)
                totalR_inSpill.append(channel_pair_spill)
                totalR_postSpill.append(channel_pair_post)
    
    
    # total energy over entire run
    congregL = getCongregateChannelData(totalL)
    congregR = getCongregateChannelData(totalR)
    
    # aggrgate energy spectra 
    L = congregL[:,1]
    single = fitAggregateSpectra(L,str(int(toGeoChannelID(mostpopulated_L))),display = True)
    
    # now use getEnergyCounts function to find total counts, counts in photopeak, and counts past photopeak (prompt gammas)
    counts = getEnergyCounts(L,single[0],single[1])
    print('Total Counts:',counts[0])
    print('Photopeak Counts:',counts[1])
    print('Prompt Counts:',counts[2])
    
    R = congregR[:,4]
    single = fitAggregateSpectra(R,int(toGeoChannelID(mostpopulated_R)),display = True)
    counts = getEnergyCounts(R,single[0],single[1])
    print('Total Counts:',counts[0])
    print('Photopeak Counts:',counts[1])
    print('Prompt Counts:',counts[2])
    
    
    # in spill energy 
    # dont fit to these
    congregL_spill = getCongregateChannelData(totalL_inSpill)
    congregR_spill = getCongregateChannelData(totalR_inSpill)
    L = congregL_spill[:,1]
    y,x,_ = plt.hist(L)
    plt.close()
    print('Total Counts: ',np.sum(y))
    R = congregR_spill[:,4]
    y,x,_ = plt.hist(R)
    plt.close()
    print('Total Counts: ',np.sum(y))
    
    
    # post-spill energy
    congregL_post = getCongregateChannelData(totalL_postSpill)
    congregR_post = getCongregateChannelData(totalR_postSpill)
    
    L = congregL_post[:,1]
    single = fitAggregateSpectra(L,int(toGeoChannelID(mostpopulated_L)),display = True)
    counts = getEnergyCounts(L,single[0],single[1])
    print('Total Counts:',counts[0])
    print('Photopeak Counts:',counts[1])
    print('Prompt Counts:',counts[2])
    
    R = congregR_post[:,4]
    single = fitAggregateSpectra(R,int(toGeoChannelID(mostpopulated_R)),display = True)
    counts = getEnergyCounts(R,single[0],single[1])
    print('Total Counts:',counts[0])
    print('Photopeak Counts:',counts[1])
    print('Prompt Counts:',counts[2])
    
