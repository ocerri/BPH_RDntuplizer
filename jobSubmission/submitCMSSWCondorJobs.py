#!/usr/bin/env python
import os, sys, subprocess, re
import argparse
import commands
import time
from glob import glob

#____________________________________________________________________________________________________________
### processing the external os commands
def processCmd(cmd, quite = 0):
    status, output = commands.getstatusoutput(cmd)
    if (status !=0 and not quite):
        print 'Error in processing command:\n   ['+cmd+']'
        print 'Output:\n   ['+output+'] \n'
    return output

#_____________________________________________________________________________________________________________
#example line: python jobSubmission/submitCMSSWCondorJobs.py -i /eos/user/o/ocerri/BPhysics/data/cmsMC_private/BPH_Tag-B0_MuNuDmst-pD0bar-kp_13TeV-pythia8_SoftQCD_PTFilter5_0p0-evtgen_HQET2_central_PU35_10-2-3_v0/jobs_out/*MINIAODSIM*.root -o /afs/cern.ch/user/o/ocerri/cernbox/BPhysics/data/cmsMC_private/BPH_Tag-B0_MuNuDmst-pD0bar-kp_13TeV-pythia8_SoftQCD_PTFilter5_0p0-evtgen_HQET2_central_PU35_10-2-3_v0/MuDst_candidates/out.root -f -c config/cmssw_privateMC_Tag_MuDmst-pD0bar-kp.py --maxtime 8h -N 5 --name MC
if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument ('-N', '--nFilePerJob', type=int, help='Number of files per job', default=10)
    parser.add_argument ('--nMaxJobs', type=int, help='Maximum number of jobs that will be submitted', default=None)
    parser.add_argument ('-i', '--input_file', type=str, default='', help='Input file template for glob or list in a txt format', nargs='+')
    parser.add_argument ('-o', '--output_file', type=str, default='', help='Output root file template')
    parser.add_argument ('-f', '--force_production', action='store_true', default=False, help='Proceed even if output file is already existing')
    parser.add_argument ('-c', '--config', type=str, help='Config file for cmsRUn')
    parser.add_argument ('--maxtime', help='Max wall run time [s=seconds, m=minutes, h=hours, d=days]', default='8h')
    parser.add_argument ('--memory', help='min virtual memory in MB', default='2000')
    parser.add_argument ('--disk', help='min disk space in MB', default='4000')
    parser.add_argument ('--cpu', help='cpu threads', default='1')
    parser.add_argument ('--name', type=str, default='BPH_RDntuplizer', help='Job batch name')

    args = parser.parse_args()

    if not sys.argv[0].startswith('jobSubmission'):
        print 'You have to run from the BPH_RDntuplizer directory'
        exit()

    '''
    ######################## Prepare output ###################################
    '''
    if not args.output_file:
        print 'No output file provided'
        exit()
    if not args.output_file.endswith('.root'):
        print 'Outputfile must end with .root'
        exit()
    outdir = os.path.dirname(args.output_file)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    else:
        if os.listdir(outdir):
            if args.force_production:
                os.system('rm -rf '+outdir+'/*')
            else:
                print 'Output directory not empty'
                print 'Empty the given directory, change directory or run with -f'
                exit()

    os.makedirs(outdir+'/out/')
    os.makedirs(outdir+'/cfg/')


    '''
    ################### Prepare input file and division #######################
    '''
    if not args.input_file:
        print 'No input file(s) provided'
        exit()

    flist = []
    if len(args.input_file) == 1 and args.input_file[0].endswith('.txt'):
        with open(args.input_file[0]) as f:
            for l in f.readlines():
                flist.append(l[:-1])
    elif len(args.input_file) == 1:
        flist = glob(args.input_file[0])
    elif isinstance(args.input_file, list):
        flist = args.input_file

    rem = len(flist)%args.nFilePerJob
    Njobs = len(flist)/args.nFilePerJob
    if rem: Njobs += 1

    print 'Input file provided:', len(flist)
    print 'Will be divided into', Njobs, 'jobs'
    for i in range(Njobs):
        i_start = i*args.nFilePerJob
        i_end = min((i+1)*args.nFilePerJob, len(flist))
        aux = '\n'.join(flist[i_start:i_end])
        with open(outdir + '/cfg/file_list_{}.txt'.format(i), 'w') as f:
            f.write(aux+'\n')

    if args.nMaxJobs:
        print 'Max number of jobs set to ', args.nMaxJobs
        if args.nMaxJobs < Njobs:
            Njobs = args.nMaxJobs
            print 'Only', Njobs, 'will be submitted'
    if Njobs == 0: exit()
    '''
    ###################### Check CMSSW and config ############################
    '''

    if not args.config:
        print 'No config file provided'
        exit()
    if not os.path.exists(args.config):
        print 'Config non existing'
        exit()

    ntuplizer_loc = os.environ['PWD']


    ''' ###################### Creating the sub #############################'''
    time_scale = {'s':1, 'm':60, 'h':60*60, 'd':60*60*24}
    maxRunTime = int(args.maxtime[:-1]) * time_scale[args.maxtime[-1]]

    os.system('chmod +x jobSubmission/CMSSWCondorJob.sh')
    print 'Creating submission scripts'

    with open('jobs.sub', 'w') as fsub:
        fsub.write('executable    = {}/jobSubmission/CMSSWCondorJob.sh\n'.format(os.environ['PWD']))

        exec_args = os.environ['PWD']+' '+args.config
        exec_args += ' ' + outdir + '/cfg/file_list_$(ProcId).txt'
        exec_args += ' ' + args.output_file.replace('.root', '_$(ProcId).root')
        fsub.write('arguments      = ' + exec_args)
        fsub.write('\n')
        fsub.write('output         = {}/out/job_$(ProcId)_$(ClusterId).out'.format(outdir))
        fsub.write('\n')
        fsub.write('error          = {}/out/job_$(ProcId)_$(ClusterId).err'.format(outdir))
        fsub.write('\n')
        fsub.write('log            = {}/out/job_$(ProcId)_$(ClusterId).log'.format(outdir))
        fsub.write('\n')
        fsub.write('+MaxRuntime    = '+str(maxRunTime))
        fsub.write('\n')
        if os.uname()[1] == 'login-1.hep.caltech.edu':
            fsub.write('+RunAsOwner = True')
            fsub.write('\n')
            fsub.write('+InteractiveUser = True')
            fsub.write('\n')
            fsub.write('+SingularityImage = "/cvmfs/singularity.opensciencegrid.org/bbockelm/cms:rhel7"')
            fsub.write('\n')
            fsub.write('+SingularityBindCVMFS = True')
            fsub.write('\n')
            fsub.write('run_as_owner = True')
            fsub.write('\n')
            fsub.write('RequestDisk = ' + args.disk)
            fsub.write('\n')
            # fsub.write('RequestMemory = ' + args.memory) #Static allocation
            fsub.write('request_memory = ifthenelse(MemoryUsage =!= undefined, MAX({{MemoryUsage + 1024, {0}}}), {0})'.format(args.memory)) # Dynamic allocation
            fsub.write('\n')
            fsub.write('RequestCpus = ' + args.cpu)
            fsub.write('\n')
        fsub.write('+JobBatchName  = '+args.name)
        fsub.write('\n')
        fsub.write('x509userproxy  = $ENV(X509_USER_PROXY)')
        fsub.write('\n')
        fsub.write('on_exit_remove = (ExitBySignal == False) && (ExitCode == 0)')
        fsub.write('\n')
        fsub.write('on_exit_hold = (ExitBySignal == True) || (ExitCode != 0)')   # Send the job to Held state on failure.
        fsub.write('\n')
        fsub.write('periodic_release =  (NumJobStarts < 3) && ((CurrentTime - EnteredCurrentStatus) > (60*60))')   # Periodically retry the jobs for 3 times with an interval 1 hour.
        fsub.write('\n')
        fsub.write('max_retries    = 3')
        fsub.write('\n')
        fsub.write('requirements   = Machine =!= LastRemoteHost')
        fsub.write('\n')
        fsub.write('universe = vanilla')
        fsub.write('\n')
        fsub.write('queue '+str(Njobs))
        fsub.write('\n')

    print 'Submitting jobs'
    output = processCmd('condor_submit jobs.sub')
    print 'Job submitted'
    os.system('mv jobs.sub '+outdir+'/cfg/jobs.sub')
