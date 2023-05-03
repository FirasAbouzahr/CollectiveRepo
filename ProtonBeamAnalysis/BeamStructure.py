from ProtonBeamHeader import *

# these lists are unique from our trip to the houston proton center on 03/05/23, so change accordingly
dir = "Flash_Therapy/PET_3-5-23/"

runLengths = [180,180,180,1200,1200,900,900,900,300,1200,600,600,900]
spillDurations = [101.7,101.5,101.6,102.0,102.0,101.7,101.6,101.8,101.8,101.6,101.6,101.8,101.3] # in ms
runNums = np.arange(1,14,1) # make sure getPaths() grabs data in order of run number
datafiles = getPaths(dir)[0]
datafiles.remove('.DS_Store')

PETx_spill_list, PETy_spill_list, PETx_list, PETy_list = coincidenceHistograms(dir,datafiles,5,5,runLengths)

postList = []
spillList = []
spillTime = [] # useful to keep when the beam turns off for energy analysis and PET imaging

# manually select which data Runs you want to look at, will plot in-beam data and post-beam
selection = [5,7,13]

for runNum in runNums:
    
    if runNum in selection:

        # isolated spill data by 250 ms on each side done via the modified Flash_DQ
        PETy = PETy_spill_list[runNum - 1]
        PETx = PETx_spill_list[runNum - 1]
        width = (max(PETx)-min(PETx))/len(PETx)

        # plot it for reference
        fig,ax = plt.subplots()
        plt.fill_between(PETx[:len(PETy)],PETy,step='mid')
        fig.set_size_inches(16,9)
        plt.show()

        #######################################################################
        #######################################################################
        ######################## Grab PET Response Data #######################
        ########################    during Beam Spill   #######################
        #######################################################################
        #######################################################################

        # look at reference plot and by eye make the pre-flash cut
        # found this method is better than an algorithim based on deviation from mean since extinction coefficient is non-zero and varies seemingly randomly
        cut = float(input('What is the pre-flash cut: '))

        for i,j in zip(range(len(PETy)),range(len(PETx))):
            if PETy[i] > cut:
                start = i
                break

        PETspilly = PETy[start:]
        PETspillx = PETx[start:]

        # grab known spill duration
        dur = spillDurations[runNum - 1]

        # temporary t, this will include time from spill to isotope decay
        t = np.linspace(0,len(PETspillx)*5,len(PETspillx))
        PETspill_t = [] # Beam Spill time
        PETspill_y = [] # Beam Spill Coincidence Counts

        # cut on known spill lengths
        for i,j in zip(t,PETspilly):
            if i > dur:
                break
            PETspill_t.append(i) # fully isolated spill lists
            PETspill_y.append(j)
        
        spillList.append([PETspill_t,PETspill_y])
        
        #######################################################################
        #######################################################################
        ######################  Grab Post-Spill PET Data  #####################
        ###################### (Activated Isotope Decays) #####################
        #######################################################################
        #######################################################################

        postx = PETx_list[runNum - 1]
        posty = PETy_list[runNum - 1]

        # grab start of post decay based on where the beam spill ended
        time_of_post = np.where(PETspilly == PETspill_y[len(PETspill_y)-1])
        time_of_post = time_of_post[0][0]

        # find where the beam end index is located in our run data
        post_start = np.where(postx <= PETspillx[time_of_post])
        post_start = post_start[0][0]
        
        # finally cut on this time to isolate post beam spill data ONLY
        post_t = postx[post_start:]
        post_y = posty[post_start:]
        
        # grab spill times, this is useful for energy response analysis
        # appending [spill start, spill end/post start]
        spillTime.append([PETspillx[0],PETspillx[0]+dur/1000])
        
        post_t = np.linspace(0,len(post_t)*.5,len(post_t))
        postList.append([post_t[:len(post_y)],post_y])
    

# now spillList & postList have completely isolated data for each region respectively, example of their plots:

fig,ax = plt.subplots(figsize = (16,9))
plt.step(spillList[0][0][1:],spillList[0][1][1:],color = 'blue',label = 'Run ' + str(runs[0]))
plt.step(spillList[1][0][1:],spillList[1][1][1:],color = 'green',label = 'Run ' + str(runs[1]))
plt.step(spillList[2][0][1:],spillList[2][1][1:],color = 'red',label = 'Run ' + str(runs[2]))
ylims = ax.get_ylim()
plt.yscale('log')
lims = list(ax.get_ylim())
plt.xticks(np.arange(10,110,10))
plt.yticks([5*10**3,7.5*10**3,10**4,2.5*10**4])
plt.legend()
plt.xlabel('Time [ms]')
plt.ylabel('No. Coincidences per 5 ms')
# plt.yscale('log')

fig,ax = plt.subplots(figsize = (16,9))
plt.step(postList[0][0][1:],postList[0][1][1:],color = 'blue',label = 'Run ' + str(runs[0]))
plt.step(postList[1][0][1:],postList[1][1][1:],color = 'green',label = 'Run ' + str(runs[1]))
plt.step(postList[2][0][1:],postList[2][1][1:],color = 'red',label = 'Run ' + str(runs[2]))
ylims = list(ax.get_ylim())
xlims = list(ax.get_xlim())
plt.yscale('log')
plt.xticks(np.arange(0,130,10))
plt.legend()
plt.xlabel('Time [s]')
plt.ylabel('No. Coincidences per 5 s')
