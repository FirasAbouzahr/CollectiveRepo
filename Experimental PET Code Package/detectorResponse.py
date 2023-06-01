from PETheader import *

# for obtaining energy/timing resolutions
def gaussian(x,A,mu,sig):
    return A * np.exp(-((x-mu)/sig)**2)

def SingleChannelEnergyResponse(df,channelID,bins):
    if channelID in np.unique(df.ChannelIDL):
        df = df[df.ChannelIDL == channelID]
        data = df.ChargeL
    else:
        df = df[df.ChannelIDR == channelID]
        data = df.ChargeR
    
    fig,ax = plt.subplots(figsize = (10,7))
    y,x,_ = plt.hist(data,bins = bins)
    bincenters = np.array([0.5 * (x[i] + x[i+1]) for i in range(len(x) - 1)])
    
    A = np.max(y)
    mu = x[np.where(y == A)[0][0]]
    guess = [A,mu,1]
    
    try:
        p,c = curve_fit(gaussian,bincenters,y,p0=guess)
        xspace = np.linspace(p[1]-2.5*p[2],p[1]+2.5*p[2])
        ax.plot(xspace,gaussian(xspace,*p),color = 'red',linestyle='dashed')
        FWHM = abs(2.355 * p[2])
        eres = FWHM / p[1] * 100
    except:
        print('Fit-Failed')
        eres = None
    return eres


def getCoincidenceTimeDiffs(df,IDL,IDR,bins):
    df_coinc = df[df.ChannelIDL == IDL]
    df_coinc = df_coinc[df_coinc.ChannelIDR == IDR]
    
    timeDiffs = []
    for timeL,timeR in zip(df_coinc.TimeL,df_coinc.TimeR):
        timeDiffs.append(timeL - timeR)
    
    fig,ax = plt.subplots(figsize = (10,7))
    y,x,_ = plt.hist(timeDiffs,bins = bins)
    bincenters = np.array([0.5 * (x[i] + x[i+1]) for i in range(len(x) - 1)])
    
    A = np.max(y)
    mu = x[np.where(y == A)[0][0]]
    std_cut = x[x >= mu - 2*mu]
    std = np.std(std_cut[std_cut <= 2*mu])
    guess = [A,mu,std]
    
    try:
        p,c = curve_fit(gaussian,bincenters,y,p0=guess)
        xspace = np.linspace(p[1]-2.5*p[2],p[1]+2.5*p[2])
        ax.plot(xspace,gaussian(xspace,*p),color = 'red',linestyle='dashed')
        CTR = abs(2.355 * p[2])
    except:
        print('Fit-Failed')
        eres = None
    return CTR
    
