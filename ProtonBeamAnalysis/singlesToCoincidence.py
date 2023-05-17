from ProtonBeamHeader.py import *

# although our readout electronics can already process data into singles, this serves as a way to do it ourselves to understand coincidence structure during beam spills better.

def SinglesToCoincidence(dfarray,timewindow):
    dfcoincidence = {'TimeL':[],'ChargeL':[],'AbsChanIDL':[],'TimeR':[],'ChargeR':[],'AbsChanIDR':[],'TimeDiff':[]}

    for hit1 in range(len(dfarray)):
        t1 = dfarray[hit1][0]
        c1 = dfarray[hit1][1]
        p1 = dfarray[hit1][2]

        if p1 in L:
            for hit2 in range(hit1+1,len(dfarray)):
                t2 = dfarray[hit2][0]
                c2 = dfarray[hit2][1]
                p2 = dfarray[hit2][2]

                if p2 in R and abs(t2-t1) <= timewindow:
                    print('Time window of',timewindow,'met. Time Diff =',abs(t2-t1))
                    dfcoincidence['TimeL'] += [t1]
                    dfcoincidence['ChargeL'] += [c1]
                    dfcoincidence['AbsChanIDL'] += [p1]
                    dfcoincidence['TimeR'] += [t2]
                    dfcoincidence['ChargeR'] += [c2]
                    dfcoincidence['AbsChanIDR'] += [p2]
                    dfcoincidence['TimeDiff'] += [abs(t2-t1)]

                else:
                    break
    
    return pd.DataFrame(dfcoincidence)

# example of usage below:

# read spill isolated singles data in

df = pd.read_csv('Desktop/COINC/HDPE1_Singles_Spill_new.dat')
df.columns = ['Time','Charge','AbsChanID']
df['Time'] = df['Time'] / 1000000000000
geoID = []
for ID in df['AbsChanID'].to_numpy():
    geoID.append(toGeoChannelID(ID))
df['AbsChanID'] = np.array(geoID).astype(int)

geo = np.split(np.sort(np.unique(geoID)),2)
L = geo[0]
R = geo[1]

dfarray = df.to_numpy()


# create custom made coincidence files for a variety of time windows
timewindow = 1e-9
wins = [1,3,5,8,10]

for i in wins:
    coincidence_df = SinglesToCoincidence(dfarray,i*1e-9)
    coincidence_df.to_csv('Coincidence_' +str(i)+'ns.txt')


# reading in the files we read out above:
df1ns = pd.read_csv('Desktop/COINC/Coincidence_1ns.txt')
df1ns.columns = ['index','TimeL','ChargeL','AbsChanIDL','TimeR','ChargeR','AbsChanIDR','TimeDiff']

df3ns = pd.read_csv('Desktop/COINC/Coincidence_3ns.txt')
df3ns.columns = ['index','TimeL','ChargeL','AbsChanIDL','TimeR','ChargeR','AbsChanIDR','TimeDiff']

df5ns = pd.read_csv('Desktop/COINC/Coincidence_5ns.txt')
df5ns.columns = ['index','TimeL','ChargeL','AbsChanIDL','TimeR','ChargeR','AbsChanIDR','TimeDiff']

df8ns = pd.read_csv('Desktop/COINC/Coincidence_8ns.txt')
df8ns.columns = ['index','TimeL','ChargeL','AbsChanIDL','TimeR','ChargeR','AbsChanIDR','TimeDiff']

df10ns = pd.read_csv('Desktop/COINC/Coincidence_10ns.txt')
df10ns.columns = ['index','TimeL','ChargeL','AbsChanIDL','TimeR','ChargeR','AbsChanIDR','TimeDiff']

# looking at the coincidence data vs time
fig,ax = plt.subplots(figsize = (16,9))
ysingles,xsingles,_ = plt.hist(df.Time,bins = 100,alpha = 0.2,label = 'Singles Data')
y1ns,x1ns,_ = plt.hist(df1ns.TimeL,bins = 100,alpha = 0.3,label = '1 ns coinc. window')
y3ns,x3ns,_ = plt.hist(df3ns.TimeL,bins = 100,alpha = 0.4,label = '3 ns coinc. window')
y5ns,x5ns,_ = plt.hist(df5ns.TimeL,bins = 100,alpha = 0.5,label = '5 ns coinc. window')
y8ns,x8ns,_ = plt.hist(df8ns.TimeL,bins = 100,alpha = 0.6,label = '8 ns coinc. window')
y10ns,x10ns,_ = plt.hist(df10ns.TimeL,bins = 100,alpha = 0.7,label = '10 ns coinc. window')
plt.yscale('log')
fig.legend()

# looking at the same coincidences vs time data but as step functions (easier to see)
fig,ax = plt.subplots(figsize = (16,9))
plt.step(xsingles[:len(ysingles)],ysingles,label = 'Singles Data')
plt.step(x1ns[:len(y1ns)],y1ns,label = '1 ns coinc. window')
plt.step(x3ns[:len(y3ns)],y3ns,label = '3 ns coinc. window')
plt.step(x5ns[:len(y5ns)],y5ns,label = '5 ns coinc. window')
plt.step(x8ns[:len(y8ns)],y8ns,label = '8 ns coinc. window')
plt.step(x10ns[:len(y10ns)],y10ns,label = '10 ns coinc. window')
# plt.yscale('log')
fig.legend(fontsize = 14,bbox_to_anchor = (1.08,.9))
plt.xticks(fontsize = 28)
plt.yticks(fontsize = 28)
plt.ylabel('')
