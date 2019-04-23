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
#example line: python submitCondorJobs.py --nev 30000 --njobs 500 --maxtime 12h --PU 0
if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument ('-N', '--nFilePerJob', type=int, help='Number of files per job', default=10)
    parser.add_argument ('-i', '--input_file', type=str, default='', help='Input file template for glob or list in a txt format', nargs='+')
    parser.add_argument ('-o', '--output_file', type=str, default='', help='Output root file template')
    parser.add_argument ('-f', '--force_production', action='store_true', default=False, help='Proceed even if output file is already existing')
    parser.add_argument ('-c', '--config', type=str, help='Config file for cmsRUn')
    parser.add_argument ('--maxtime', help='Max wall run time [s=seconds, m=minutes, h=hours, d=days]', default='8h')

    args = parser.parse_args()

    if not sys.argv[0].startswith('jobSubmission'):
        print 'You have to run from the BPH_RDntuplizer directory'
        exit()


    '''
    ################### Prepare input file and division #######################
    '''
    if not args.input_file:
        print 'No input file provided'
        exit()
    flist = []
    if isinstance(args.input_file, basestring) and args.input_file.endswith('.txt'):
        with open(args.input_file) as f:
            for l in f.readlines():
                flist.append(l)
    elif isinstance(args.input_file, basestring):
        flist = glob(args.input_file)
    elif isinstance(args.input_file, list):
        flist = args.input_file

    rem = len(flist)%args.nFilePerJob
    Njobs = len(flist)/args.nFilePerJob
    if rem: Njobs += 1

    infilestr = []
    for i in range(Njobs):
        i_start = i*args.nFilePerJob
        i_end = min((i+1)*args.nFilePerJob, len(flist))
        aux = ', '.join(flist[i_start:i_end])
        # print aux
        infilestr.append(aux)




    '''
    ######################## Prepare input ####################################
    '''
    if not args.output_file:
        print 'No output file provided'
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
                print 'Empty the give directory, change directory or run with -f'
                exit()

    os.makedirs(outdir+'/out/')
    os.makedirs(outdir+'/cfg/')

    outfiletemplate = os.path.basename(args.output_file)


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

    os.system('chmod +x jobSubmission/CMSSECondorJob.sh')
    print 'Creating submission script\n\n'

    fsub = open('jobs.sub', 'w')
    exec_base = 'executable    = {}/jobSubmission/CMSSECondorJob.sh\n'.format(os.environ['PWD'])
    fsub.write(exec_base+'job1023_gen_wPU_v1.sh\n')
    # exec_args = str(nev)+' '+str(st_seed)+' $(ProcId) '+args.process+' '+version+' '+args.CMSSW_loc
    # if isinstance(args.PU, int) and args.PU > 0:
    #     exec_args += ' ' + str(args.PU)
    # fsub.write('arguments     = ' + exec_args)
    # fsub.write('\n')
    # fsub.write('output        = {}/out/{}.$(ClusterId).$(ProcId).out'.format(outdir, args.process))
    # fsub.write('\n')
    # fsub.write('error         = {}/out/{}.$(ClusterId).$(ProcId).err'.format(outdir, args.process))
    # fsub.write('\n')
    # fsub.write('log           = {}/out/{}.$(ClusterId).$(ProcId).log'.format(outdir, args.process))
    # fsub.write('\n')
    # fsub.write('+MaxRuntime   = '+str(maxRunTime))
    # fsub.write('\n')
    # fsub.write('+JobBatchName = '+args.process)
    # fsub.write('\n')
    # fsub.write('x509userproxy = $ENV(X509_USER_PROXY)')
    # fsub.write('\n')
    # fsub.write('queue '+str(njobs))
    # fsub.write('\n')
    # fsub.close()
    #
    # output = processCmd('condor_submit jobs.sub')
    # os.rename('jobs.sub', outdir+'/cfg/jobs.sub')
