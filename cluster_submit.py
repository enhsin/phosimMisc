#!/usr/bin/env python
##
## @package phosim.py
## @file cluster_submit.py
## @brief script to submit phosim jobs to a cluster
##
## @brief Created by:
## @author Glenn Sembroski (Purdue)
##
## @brief Modified by:
##
## @warning This code is not fully validated
## and not ready for full release.  Please
## treat results with caution.
##
## The intention of this script is to use the dagMan and  .submit files
## created in the phosim "work" directory by the running of phosim with the 
##  -g condor option, to generate and submit the slurm (or pbs) files to run 
## the trim, raytrace and e2adc jobs on a cluster. And then to submit them.
## Note: Particular clusters required specific modifications to the submission
## files. (Provisions and examples of this exist in functions
## setupForHost() and getPid(0 functions).
#

import os
import sys, optparse
import subprocess

## print the usage
def usage():
     script=os.path.abspath(__file__)
     os.system(script+' x --help')

## Seach the dagMan file for the dependancy lines. They start with the word
# "PARENT". Create and return a dict of Child:Parent for use later is
# deriveing dependancy commands for the submition file
def getDependancies(dagManFileFull):
     #search the dagMan file for the "PARENT" lines. On these lines are the
     #job dependancies

     #Get all the lines
     dFile = open(dagManFileFull,"r")
     dagMan = dFile.readlines()
     dFile.close()

     dependancies = {}
     for line in dagMan:
          lstr1 = line.split()
          #see if this line is a depancncy line
          if lstr1[0] == "PARENT" and lstr1[2] == "CHILD":
               #extract things from this line
               parentJobName = lstr1[1]
               childJobName = lstr1[3]
               dependancies[childJobName] = parentJobName
     return dependancies

##getJobs returns a dict with the jobName as key and submission file name
# as the value
def getJobs(dagManFileFull,jobType):
    #Go through the dagman file and find all "JOB" lines with jobType
    #Get all the lines (I know this duplicates the reading of the file we did
    #in getDependancies but I wanted to keep thing independant and the code
    #obvious
     
    dFile = open(dagManFileFull,"r")
    dagMan = dFile.readlines()
    dFile.close()
     
    # search for "JOB" and jobType lines
    jobs = {}
    jobTypeLen=len(jobType)
    for line in dagMan:
        lstr1 = line.split()
        #see if this is a JOB and jobType line
        if lstr1[0] == "JOB":
            jobName= lstr1[1]
            if jobName[:jobTypeLen] == jobType:
                submitFileName = lstr1[2]
	        jobs[jobName] =  submitFileName
     
    return jobs

##Read througth the submision file for this job and extract the values into
# a dict.
def getSubmissionParams(submitFileName):

    sFile = open(submitFileName,"r")
    jobSubmitLines = sFile.readlines()
    sFile.close()
               
    jobSubDict = {}
    for line in jobSubmitLines:
        submitLine=line.split()
        if submitLine[0] == "executable":
            jobSubDict['executable'] = submitLine[2]
        elif submitLine[0] == "initialdir":
            jobSubDict['initialDir'] = submitLine[2]
        elif submitLine[0] == "Input":
            jobSubDict['input'] = submitLine[2]
        elif submitLine[0] == "Output":
            jobSubDict['output'] = submitLine[2]
        elif submitLine[0] == "Log":
            jobSubDict['log'] = submitLine[2]
        elif submitLine[0] == "Error":
            jobSubDict['error'] =  submitLine[2]
        elif submitLine[0] == "transfer_input_files":
            jobTransferFiles = " "
            jobTransferFiles = jobTransferFiles.join(submitLine[2:])
            jobSubDict['transferFiles'] = jobTransferFiles

    return jobSubDict

## Setup a string of pbs/slrumn values;
def setupForHost():
    host = os.getenv('HOSTNAME')
    #Setup the commands for the cluster we are using:
    #Setup the nodes and cpu per node. We will change this a little when we 
    #add the possibility of threads. 2 possible options either seperate nodes 
    #and ppn or combined(NERSC)
    #   
      
    #For NERSC we may want to list the disk we are using.
    #Returend is: base of PBS/slrumn submit command,Dependency command,
    #string of all other commands, string of maximum thread count.
    if host[:6] == 'edison':
        pbsSetupList = ('#SBATCH -L SCRATCH' + "\n" +
                        '#SBATCH -p shared' + "\n" +
                        '#SBATCH -t 10:00:00' + "\n" +
                        '#SBATCH -N 1' + "\n" +
                        '#SBATCH --mem=5GB' + "\n" +
                        '#SBATCH -A m1727'+ "\n")

        return {'SUBMITCMD':   'sbatch',
                'DEPENDCMD':   '#SBATCH -d afterok:', 
                'INITIALLIST': pbsSetupList,
	        'MAXTHREADS':  '24',
                'THREADCMD':   '#SBATCH -n '}

    elif host[:4] == 'cori':
        pbsSetupList = ('#SBATCH -L SCRATCH' + "\n" +
                        '#SBATCH -p shared' + "\n" +
                        '#SBATCH -t 10:00:00' + "\n" +
                        '#SBATCH -N 1' + "\n" +
                        '#SBATCH -C haswell' + "\n" +
                        '#SBATCH -A m1727'+ "\n")
 
        return {'SUBMITCMD':   'sbatch',
                'DEPENDCMD':   '#SBATCH -d afterok:', 
                'INITIALLIST': pbsSetupList,
                'MAXTHREADS':  "32",
                'THREADCMD':   '#SBATCH -n '}

    elif host[:6] == 'hammer' or host[:5] == 'conte':
        pbsSetupList = ('#PBS -q standby' + "\n" +
                        '#PBS -l walltime=4:00:00' + "\n" +
                        '#PBS -l mem=5GB' + "\n" +
                        '#PBS -l naccesspolicy=singleuser'+ "\n")

        return {'SUBMITCMD':   'qsub -V ',
                'DEPENDCMD':   '#PBS -W depend=afterok:', 
                'INITIALLIST': pbsSetupList,'MAXTHREADS':"20",
                'THREADCMD':   '#PBS -l nodes=1:ppn='}

    else:
        print('Unknow host: ' + host )
        print('Known hosts are: NERSC:Edison,Cori (haswell), Purdue:Hammer')
        sys.exit(1)

##Search Input pars file for "thread" command and get the number of requested
#threads (as a string). If not there return requested number of "1"
def getNumRequestedThreads(inputFileName):
    iFile = open(inputFileName,'r')
    iLines = iFile.readlines()
    iFile.close()
    thread=1
    for line in iLines:
        lstr1 = line.split()
        #see if this line is a "thread" line
        if lstr1[0] == "thread" :
            #extract thread count  from this line
            threadCount = lstr1[1]
            break

    return threadCount



## Get the cluster assiged PID for this submitted job.
#This usually only called for the trim jobns and used when setting up 
#dependencies for the raytrace /e2adc jobs.
def getPid(jobPBS,subPBSLine):
    print ('Submission of ' + jobPBS  + ' returned:' + 
           subPBSLine)
    host = os.getenv('HOSTNAME')
    if host[:6] == 'edison':
        #return from the submission will be something like:
        # "Submitted batch job 2162058". So we want the 4th word

        JobPid = subPBSLine.split()[3]

    if host[:4] == 'cori':
        #return from the submission will be something like:
        # "Submitted batch job 2162058". So we want the 4th word

        JobPid = subPBSLine.split()[3]

    if host[:6] == 'hammer' or host[:5] == 'conte':
        #return from the submission will be something like:
        # "1029466.hammer-adm.rcac.purdue.edu". So we want the string before 
        #the "."

        JobPid = subPBSLine.split(".")[0]
   
    return JobPid

## Add the sbatch (or pbs) cluster commands to the submission file.
def addClusterCommands(pfile, submitPBSList, threadCMD, numThreads):
    #Print out all the PBS/slrumn commands in the List (except for the 
    #DEPEND cmd. Thats later.
    #The submitPBSList and threadCMD strings come from the call to 
    #setupForHost 
    #Setup start of pbs submission file

    pfile.write("#!/bin/bash -l\n")
    pfile.write(submitPBSList)
    pfile.write(threadCMD + numThreads + "\n")
    return

## Determine if this is "CHILD" job. That is if it depends on a parent job
# running first. If so add a depenacny line to the submission file.
def addDependancyCommand(pfile, jobName, dependancies, jobIDs, 
                         submitDependCmd):
    # First check to see that this jobName is a child to someones parent
    # Then check that the parent has already been run and we have a pid for it
    #
    # This could be more general by not assuming the above ordering by
    # searching for jobs with no child tasks, running them first and so on
    # and so on. That would complicate thing highly but its doable
    if dependancies.get(jobName) is not None:
        parentJobName=dependancies[jobName]           #Parent trim job name

        #Now make sure that this parent has run and if so get its pid
        if jobIDs.get(parentJobName) is not None:

            pid =  jobIDs[parentJobName]
            # Add the dependancy command
            pfile.write(submitDependCmd + pid + '\n')
        else:
            print('No pid avaiable for dependancy ' + parentJobName +
                  ' needed by ' +jobName )
            sys.exit(1)
    return

## Add to the pbs file the lines to run this task. (Originally from the dagMan 
# file.
def addTaskRunLine(pfile,jobSubDict,taskName):
    #run the job with its arguments. This will be multiline with 
    #continuations.
    #The log file will need to manipulated  a bit to start with the taskname
    path,logFileName = os.path.split(jobSubDict["log"])
    logFileFull = (jobSubDict["initialDir"] + "/" + path + '/' + taskName + logFileName[3:])
    pfile.write(jobSubDict['executable'] + " \\"  + '\n')
    pfile.write("< " + jobSubDict["input"] + " \\"  + '\n')
    pfile.write("> " + logFileFull + '\n')

    return

##Add a line to the .pbs script file to remove a file. 
def addRemoveFileLine(pfile,fileName): 
    pfile.write('rm ' + fileName + '\n')
    return

##Add the commands to the submission job script that runs the job.
# Not sure about the transfer stuff for now. Wait and see whats missing
def addJobCommands(pfile,jobSubDict,taskName,jobName):
    #Add the bash commands to run this job.
    #First cd to the initial directory
    pfile.write('date' + '\n')
    pfile.write('cd ' + jobSubDict['initialDir'] + '\n')

    #If this is a raytrace file we first need to append the trimcatalog to the
    #Raytrace*.pars
    if taskName == 'raytrace' and '_0' in jobName:
        #we are going to put in a line to append the trimcatalog (which 
        #the trim jobs make to the end of the raytrace.pars file
        # Trim catalog and raytrace.pars files will be in local directory(work)
        trimCatalogFileName = 'trimcatalog_' + jobName[9:-7] + '.pars'
        raytraceParsFileName  = jobSubDict['input']
        pfile.write( 'cat ' + trimCatalogFileName + ' >> ' +
	             raytraceParsFileName + '\n')

    #I may need to transfer or make links to some files here. Not sure
    addTaskRunLine(pfile,jobSubDict,taskName)
    return


## removeFile deletes files (if they do not exist, it will catch the OSError 
# exception and silently proceed.)
def removeFile(filename):
    try:
         os.remove(filename)
    except OSError:
         pass

## Add a line to a .pbs file to mv a file to a directory
def addMoveFileLine(pfile,fileName,destDir):
    fid = fileName[7:-13]
    eid = fileName[-12:-8]
    tarName = 'lsst_'+fid+'_'+eid+'.tar'
    pfile.write('tar cf '+tarName+' *'+fid+'*'+eid+'.fits.gz --remove-files\n')
    pfile.write( 'mv ' + tarName+ ' ' + destDir + '/\n')
    #pfile.write( 'mv ' + fileName + ' ' + destDir + '/\n')
   
## Seach the dagMan file for the jobs. Create the slurm ( or .pbs
# depending on cluster, assuming NERSC edison for now)
# submit the jobs and save JobNames and cluster Job_ID in a list which is
# returned. This list will be used for dependancies
def createAndSubmitJobs(opt,dagManFileFull,dependancies):
     
    submitPBSList = setupForHost()

    #search the dagMan file for the  jobs.
    #getJobs returns a dict with the jobName as key and submission file 
    #name as the value

    trimJobDict = getJobs(dagManFileFull,"trim_")
    raytraceJobDict = getJobs(dagManFileFull,"raytrace_")
    e2adcJobDict =  getJobs(dagManFileFull,"e2adc_")

    #trimJobID is a dict used for generating dependancies. It links JobName
    #and submitted job Pid's (process id). This is created from trim jobs 
    #submissions only and used only for setting dependiencies in the 
    #combined raytrace/e2adc jobs
    #The trim jobs will have no dependencies and each 
    #raytrace job will only depend on a particular trim job (many raytrace 
    #jobs will depend on a particulat trim job). Each  e2adc 
    #job depends on a particular raytrace. We implitly enforce that 
    #dependency by including the running of the e2adc task in its 
    #Parent raytrace pbs job file.
    trimJobID={}

    #Trim jobs first(Preserve order just for looks)
    jobCount = 0
    for jobName,trimSubmitFileName in sorted(trimJobDict.items()):	
        #Get the important lines from the submission file
	jobSubDict = getSubmissionParams( trimSubmitFileName)
              
        #We have everything we need to submit this task.
        #time to build the .pbs file
        # ##############################
        #Note For the purposes of this script the .pbs abreviation is used 
        #for BOTH pbs AND  slurm based clusters)
        # ##############################

        jobPBS = opt.workDir + "/" + jobName + '.pbs'
        pfile=open(jobPBS,'w') 
          
        addClusterCommands(pfile, submitPBSList['INITIALLIST'],
                           submitPBSList['THREADCMD'], "1")
        
        addJobCommands(pfile, jobSubDict,'log',jobName)  
        #The 'log' just retains trim log file name as listed in the .submit 
        #file.  That name alteady has 'trim' in it

        #cleanup(remove .pars file from "work/" Disabled for now.
	#addRemoveFileLine(pfile, jobSubDict['input'])
                    
        # The next line will probably work but lets wait a bit before we enable
        # it. Disabled for now.
	#addRemoveFileLine(pfile,jobPBS)

        #And thats all for trim
        pfile.close()
          
        #submit the job
        #Its important here to get the Job_ID (PID) that gets returned to
        #standard output for use later in setting the raytrace job
        #dependency.
          
        submitCmd = (submitPBSList['SUBMITCMD'] + 
                    ' -e ' + jobSubDict['initialDir'] + '/errors/' + jobName +
                    '.pbs.err ' + 
                    ' -o ' + jobSubDict['initialDir'] + '/logs/' + jobName + 
                    '.pbs.log ' + jobPBS)

        print ('Submit cmd: ' + submitCmd) 

        #This is it! RUN iT!!!
        # #################################################
        
        subPBSLine = subprocess.check_output(submitCmd, shell=True)
        
        # #################################################
               
        #Debug lines follows
        #jobCount = jobCount +1
        #subPBSLine = ('Submitted batch job ' + str(jobCount) )

        #Pick up the PID for dependencies later
        # We should probably check the return code here. Later!
        JobPid = getPid(jobPBS,subPBSLine)

        #save this in a dict with the trim job name as a key. Used
        # when setting up raytrace dependencies
        trimJobID[jobName]=JobPid

        #cleanup
        #removeFile(trimSubmitFileName) Disabled for now.
               
    #Trim jobs are all now all running and we have the dependencies.
    #Build and submit the combined raytrace/e2adc jobs.

    for jobName,raytraceSubmitFileName in sorted(raytraceJobDict.items()):
        #Get the important lines form the submission file
	jobSubDict = getSubmissionParams(raytraceSubmitFileName)

        #get the requested threads(if asked for or max if exceeded)
	inputDirFull = (jobSubDict["initialDir"] + "/" +  jobSubDict["input"])
	requestedThreads = getNumRequestedThreads( inputDirFull )
        if int(requestedThreads) > int(submitPBSList['MAXTHREADS']) :
            requestedThreads=submitPBSList['MAXTHREADS']

        jobPBS = opt.workDir + "/" + jobName + '.pbs'
        pfile=open(jobPBS,'w') #Don't forget to close this file when done
          
        addClusterCommands(pfile, submitPBSList['INITIALLIST'],
                           submitPBSList['THREADCMD'], requestedThreads)
	addDependancyCommand(pfile,jobName, dependancies, trimJobID,
	                     submitPBSList['DEPENDCMD'])
        addJobCommands(pfile, jobSubDict,'raytrace',jobName)
        #addRemoveFileLine(pfile,jobSubDict['input']) #Disabled for now
 
	
        #Now add the e2adc command. Get e2adc JobName from the raytrace
        #JobName
        e2adcJobName = 'e2adc_' + jobName[9:-2]

        #Make sure we have this key in our e2adc dict and get the
        # e2adc submission file name (We could also check to see if we had 
        # the dependency on the raytrace here in dependencies also as a 
        # sanity check.)

	if e2adcJobName not in e2adcJobDict:
            print ( 'Key for e2adc job:' +  e2adcJobName + 
                    ' not in dagMan File' )
            sys.exit(1)
        if jobName == dependancies[e2adcJobName]:
            e2adcSubmitFileName = e2adcJobDict[e2adcJobName]
            e2adcSubDict = getSubmissionParams(e2adcSubmitFileName)
            addTaskRunLine(pfile, e2adcSubDict,'e2adc') #will get executable path
            #addRemoveFileLine(pfile,e2adcSubDict['input']) #Disabled for now

	#At this point we need to move the image fits files to the output
        #directory. First the completed "lsst_e" file name from tne .submit 
        #file transfer list. (other way to get this name is to make it up from
        # the jobName and to search the e2adc*.pars file for the filter line
        # This is easier.)

    	    transferList=e2adcSubDict['transferFiles'].split()
    	    for transFile in transferList:
                if transFile[:6] == "lsst_e":
                    lsst_e_FileName = transFile
                    addMoveFileLine(pfile,lsst_e_FileName,opt.outputDir)
                    break
                #Not really necessatry yto get these _a_ files but its easy so
                #do it we can remove later if never useful. Use a wild card to
                # get them all
	        #lsst_a_FileName = ("lsst_a" +  lsst_e_FileName[6:-12] + '*' +
            #                       ".fits.gz")
	        #addMoveFileLine(pfile,lsst_a_FileName,opt.outputDir)
        # The next line will probably work but lets wait a bit before we enable
        # it.
        #addRemoveFileLine(pfile,jobPBS) #Disabled for now

        #And thats all for the raytrace/e2adc .pbs file. Add date print and
        #close it!
        pfile.write('date' + '\n')
        pfile.close()

        #now we can submit it. Don't need the pid. Nobody will depend on
        #this(However I may take this back if we want to add an end-of-run 
        #cleanup task which will have to wait for all the raytrace/e2adc jobs
        #to finish)
 
        submitCmd = ( submitPBSList['SUBMITCMD'] + 
                      ' -e ' + jobSubDict['initialDir'] + '/errors/' + 
                      jobName + '.pbs.err ' + 
                      ' -o ' + jobSubDict['initialDir']  +  '/logs/' + 
                      jobName + 'pbs.log ' +   jobPBS)
  
        print ('Submit cmd: ' + submitCmd) 

        # #################################################
        #status = subprocess.call(submitCmd, shell=True)
        subPBSLine = subprocess.check_output(submitCmd, shell=True)

        # #################################################

        JobPid = getPid(jobPBS,subPBSLine)
        trimJobID[jobName]=JobPid

               
	#Debug lines follows
        #jobCount = jobCount +1
        #subPBSLine = ('Submitted batch job ' + str(jobCount) )
                    
        print ( jobName + ': ' + subPBSLine)
        #We should probably check the status here (did command succdeed?)
        #Later!

        #cleanup
	#removeFile(raytraceSubmitFileName)
        #removeFile(e2adcSubmitFileName)

    return

##main function
# First we will only be interested in getting the trim jobs running. We will
# save the trim  submission job ID's for later use as a dependency requirment
# for the raytrace jobs
def main():

     print(sys.argv[0] + ': v1.1')

     # Get the current directories as defaults. This sssumes we are running
     # from the tools direcoty. These may all be overridden by commandl line
     # arguments
     toolsDir=os.path.split(os.path.abspath(__file__))[0]
     phosimDir = (toolsDir + "/..")
     defaultOutputDir=os.getenv('PHOSIM_OUTPUT_DIR', phosimDir+'/output')
     defaultWorkDir=os.getenv('PHOSIM_WORK_DIR', phosimDir+'/work')

     # Parse the command line
     parser = optparse.OptionParser(usage='%prog dagManFile ' +
	                            ' [<arg1> <arg2> ...]')

     # define acceptable options
     parser.add_option('-w','--work',dest="workDir",default=defaultWorkDir,
                       help='temporary work directory')
     parser.add_option('-o','--output',dest="outputDir",
                       default=defaultOutputDir, help='output directory')
     
     if len(sys.argv)<2:
          usage()
          sys.exit()
     if sys.argv[1] in ('-h', '--help'):
          usage()
          sys.exit()
               
     # parse_args returns a pair of values, remainder should have dagMan
     # file name
     opt, remainder = parser.parse_args(sys.argv[1:]) 

     #require dagMan file path/name
     if len(remainder) == 0:
          print 'dagMan File name argument is required'
          usage()
          sys.exit()
     if len(remainder) != 1:
          print ('Too many arguments on command line: %s',  remainder)
          usage()
          sys.exit()
     
     dagManFile = remainder[0]
     # see if we have a path to this file or if we should assume the work
     #directory
     dagDir=os.path.split(dagManFile)[0]
     if len(dagDir) == 0:
          dagManFileFull = opt.workDir + '/' + dagManFile
     else:
          dagManFileFull = dagManFile
     
     #Probably should Make sure paths to work/logs and work/errors

     # Entire workflow for this program follows
     #First create and submit the trim jobs. Returns a dict with trim Job name
     #as key and job_ID as value
     
     dependancies = getDependancies(dagManFileFull)
 
     createAndSubmitJobs(opt, dagManFileFull, dependancies)
          
if __name__ == "__main__":
    main()
