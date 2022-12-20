import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

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
def getData(df):
    return df.Area.to_numpy() * 10**9

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
# datalist = [array of data arrays]
# namelist = [data labels in same order they occur in datalist]
def getSpectra(datalist,namelist,title,bins):
    fig0,ax0 = plt.subplots()
    for data,names,t in zip(datalist,namelist,times):
        fig1,ax1 = plt.subplots()
        y,x,_ = ax1.hist(data,bins = bins)
        counts = str(int(np.sum(y)))
        ax0.step(x[:len(y)],y,label = names + ' (Counts: ' + counts +')')
        plt.title(title + " (" + (names) + ")",fontsize = 20)
        plt.ylabel('Counts',fontsize = 20)
        plt.xlabel("Pulse Area (nVs)",fontsize = 20)
        plt.xticks(fontsize = 13)
        plt.yticks(fontsize = 13)
        fig1.set_size_inches(10,7)
    
    ax0.legend(fontsize = 15)
    ax0.set_title(title,fontsize = 20)
    ax0.set_ylabel('Counts',fontsize = 20)
    ax0.set_xlabel('Pulse Area (nVs)',fontsize = 20)
    plt.xticks(fontsize = 13)
    plt.yticks(fontsize = 13)
    fig0.set_size_inches(10,7)
Footer
© 2022 GitHub, Inc.
Footer navigation
Terms
Privacy
Security
Status
Docs
Contact GitHub
Pricing
API
Training
Blog
About
