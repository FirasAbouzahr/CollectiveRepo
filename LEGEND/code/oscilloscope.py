import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from os import walk

# since we work with many data files at a time, this function grabs all file names from a directory 
# so we do not have to type these out by hand!
def getPaths(directory):
    pathlist = []
    for (dirpath, dirnames, filenames) in walk(directory):
        pathlist.append(filenames)
        break
    return pathlist
    
'''
The function removeHeader is used to remove the header from oscilloscope data files. 
MORE SPECIFICALLY it removes the first 4 lines of the data files.
This means that if you appply this function twice you will truncate data...
hence only use this once EVER per file.
'''
def removeHeader(path):
    with open(path, "r+") as file:
        cut = file.readlines()
        file.seek(0)
        for i in range(len(cut)):
            if i >= 5:
                file.write(cut[i])
        file.truncate()

# the oscilloscope file only come in one format (for our purposes), 
# so this function takes care of importing the data
def getFrame(path):
    df = pd.read_csv(path,header = None,usecols = [1])
    df.columns = ['Area']
    return df

# the function takes the result of getFrame and creates a usable array. 
# hence if you want data you do Area_example = getData(getFrame(/your/path/to/the/data))
def getData(df,conversion):
    return df.Area.to_numpy() * conversion

# this is very simple yet very important, it returns custom bins based on a specific set of data
def getBins(data,numbins):
    return np.linspace(data.min(),data.max(),numbins)

# normalize to a specific time
# realtime is the length of time the data was taken over
# timenorm is the time you would like to truncate it at
def timeNormalize(data,realtime,timenorm):
    rate = len(data)/realtime
    datacut = int(rate*timenorm)
    return data[:datacut]

# self explanatory... note we drop negatives.
def subtractBackground(data,bg,bins):
    yData,xData,_ = plt.hist(data,bins = bins)
    plt.close()
    yBKG,xBKG,_ = plt.hist(bg,bins = bins)
    plt.close()
    y = yData - yBKG
    y[y < 0] = 0
    bg_subtracted_data =  [xData[:len(yData)],y]
    return bg_subtracted_data
    
# this takes multiple sets of data and creates individual plots for each
# and creates one cumulative plot 
def getSpectra(datalist,namelist,bins,title,xlabel = '',ylabel = '',customlimit = False,lims = [0,1]):
    fig0,ax0 = plt.subplots()
    for data,names in zip(datalist,namelist):
        fig1,ax1 = plt.subplots()
        y,x,_ = ax1.hist(data,bins = bins)
        counts = str(int(np.sum(y)))
        avg = str(np.round(np.mean(data),4))
        ax0.step(x[:len(y)],y,label = names + ' (Counts: ' + counts +' | Mean: ' + avg + ')')
        plt.title(title + " (" + (names) + ")",fontsize = 20)
        plt.ylabel(ylabel,fontsize = 20)
        plt.xlabel(xlabel,fontsize = 20)
        plt.xticks(fontsize = 13)
        plt.yticks(fontsize = 13)
        fig1.set_size_inches(10,7)
        ax1.set_ylim(0,8000)
        if customlimit == True:
            plt.xlim(lims[0],lims[1])
    
    ax0.legend(fontsize = 12)
    ax0.set_title(title,fontsize = 20)
    ax0.set_ylabel(ylabel,fontsize = 20)
    ax0.set_xlabel(xlabel,fontsize = 20)
    plt.xticks(fontsize = 13)
    plt.yticks(fontsize = 13)
    fig0.set_size_inches(10,7) 
    ax0.set_ylim(0,8000)
    if customlimit == True:
        ax0.set_xlim(lims[0],lims[1])

def getMeans(datalist,roundoff = True,decimals = 3):
    avgs = []
    for i in datalist:
        if roundoff == True:
            avg = np.round(np.mean(i),decimals)
            avgs.append(avg)
        else:
            avgs.append(np.mean(i))
    return np.array(avgs)

def getCounts(datalist):
    countList = []
    for i in datalist:
        countList.append(int(len(i)))
    return np.array(countList)

def addlabels(x,y,fontsize):
    for i in range(len(x)):
        plt.text(i-.2,y[i],y[i],fontsize = fontsize)
        
def getStatComparison(statlist,namelist,title,xlabel = '',ylabel = ''):
    x = np.arange(0,len(statlist))
    fig,ax = plt.subplots()
    ax.bar(x,statlist)
    plt.xticks(x,namelist,fontsize = 13,rotation = 'vertical')
    plt.yticks(fontsize = 13)
    plt.title(title,fontsize = 20)
    plt.xlabel(xlabel,fontsize = 20)
    plt.ylabel(ylabel,fontsize = 20)
    fig.set_size_inches(10,7)
    addlabels(x,statlist)
