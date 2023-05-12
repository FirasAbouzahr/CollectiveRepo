# Firas Abouzahr

# do not include file type extensions for file names

def launcher_file(fileName, # custom launcher file name
                  executable, # name of simulation executable 
                  macroFile, # name of macrofile that tells Geant4 what to execute in simulation
                  outputFiles, # name of files to output data to 
                  numInParallel): # number of simulations you want to run in parallel
    
    file = open(fileName, "w")
    
    for nums in range(numInParallel):
        outputString = ''
        for outs in outputFiles:
            outputString += outs + str(nums) + '.txt '
            
        file.write('./' + executable + ' ' + macroFile + '.mac ' + outputString + '\n')

def job_submission(bashName, # shell script file name
                   jobName,  # name of the job, just for reference
                   errorName, # log file to check for status & errors during job
                   numNodes, # number of requested nodes to use 
                   numInParallel, # number of parallel jobs per node
                   queue, # respective supercomputer queue
                   time, # run time of job (hh:mm:ss)
                   allocationName, # name of allocation to charge job to 
                   launcherDirectory, # directory where launcher file is located
                   launcherName): # launcher file name
    
    file = open(bashName + '.bash', "w")
    
    commands = ['-J','-o','-N','-n','-p','-t','-A']
    params = [jobName,errorName,str(numNodes),str(numInParallel),queue,time,allocationName]
    file.write('#!/bin/bash' + '\n')
    
    for i,j in zip(commands,params):
        file.write('#SBATCH ' + i + ' ' + j + '\n')
    
    file.write('\n')
    file.write('module load launcher' + '\n')
    file.write('\n')
    file.write('export LAUNCHER_WORKDIR=' + launcherDirectory + '\n')
    file.write('export LAUNCHER_JOB_FILE=' + launcherName + '\n')
    file.write('\n')
    file.write('${LAUNCHER_DIR}/paramrun')

# example: generates files needed to run 10 PET simulations in parallel on TACC in the development queue
# see included output files in github generated from this example
launcher_file('launcher','PETSim','Dung',['electron','gamma','scintillation','detected'],5)
job_submission('serialJob','PETSimulation','launcher.o%j',1,20,'development','00:60:00','PET','/home1/08038/firas/sparse/build','launcher')
