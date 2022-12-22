import pandas as pd
import numpy as np

'''
Firas Abouzahr 

This code converts experimental spectrometer data to Geant4 Simulation commandline files allowing us to
directly use our labratory light sources in simulation.
   
'''

#######################################
##### Raw Data to Histogram Data ######
#######################################

dataframe = pd.read_csv('Desktop/VUVRun3.csv',usecols = [0,1,4,7,10,13]) # change usecols list depending on .csv
dataframe.columns = ['wavelength','hour_0','hour_2','hour_19','hour_21','hour_42'] # change column names accordingly
dataframe = dataframe.drop(labels = [0,1,2,3,4])
dataframe['wavelength'] = pd.to_numeric(dataframe['wavelength'],errors='coerce')

minWL = int(np.floor(dataframe.wavelength.min())) + 1
maxWL = int(np.floor(dataframe.wavelength.max()))

dict = {}
for labels in dataframe.columns:
    dict[labels] = []

for i in range(minWL,maxWL,1):
    wl = dataframe[dataframe['wavelength'] < i]
    wl = wl[wl['wavelength'] > i - 1]
    for j in wl.columns:
        if j == 'wavelength':
            dict['wavelength'] += [i - 1]
        else:
            total = int(wl[j].sum())
            dict[j].append(total)

histogramData = pd.DataFrame(dict)
histogramData.to_csv('HistogramData.csv', index=False)


#################################################
####### Histogram Data to G4 input command ######
#################################################

column = histogramData.columns[5] # change this to select what column count column you want to write out

def wavelength_to_Energy(wavelength):
    E = 1240/wavelength
    return E

g4Array = [np.array(histogramData.wavelength)]
for labels in histogramData.columns:
    if labels == column:
        x = np.array(histogramData[labels])
        g4Array.append(x)

g4File = open('g4input_' + column + '.txt', 'w') # add a prefix before g4_input to save file in a specific place

for wavelength,events in zip(g4Array[0],g4Array[1]):
    energy = round(wavelength_to_Energy(wavelength),3)
    g4File.write('/gps/ene/mono ' + str(energy) + ' eV\n')
    g4File.write('/run/beamOn ' + str(events) + '\n')
