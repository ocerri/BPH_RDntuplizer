#!/usr/bin/env python
import os, sys, subprocess, re
import argparse
import commands
import time
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

    parser.add_argument ('--njobs', help='Number of jobs in which split the production', default=10)
    parser.add_argument("-c", "--config", type=str, default=None, help="Config file for cmsRUn")
    parser.add_argument("-i", "--input_file", type=str, default=None, help="input root files", nargs='+')
    parser.add_argument("-o", "--output_file", type=str, default=None, help="output root file")

    parser.add_argument ('--CMSSW_loc', help='CMSSW src loc', default='/afs/cern.ch/user/o/ocerri/cernbox/BPhysics/CMSSW_10_2_3/src/')
    parser.add_argument ('--force_production', action='store_true', default=False, help='Proceed even if output file is already existing')
    parser.add_argument ('--maxtime', help='Max wall run time [s=seconds, m=minutes, h=hours, d=days]', default='8h')

    args = parser.parse_args()

    print 'Done up to here. To be ultimated.'
    print 'Only RD config is ready for this submission'
    exit()

    nev        = args.nev
    njobs      = int(args.njobs)
    st_seed    = int(args.st_seed)
    cmssw_version = re.search('/CMSSW_[0-9]+_[0-9]+_[0-9]+/src', args.CMSSW_loc).group(0)[7:-4]
    version    = cmssw_version.replace('_', '-') + '_' + args.version
    if args.PU == 0:
        version = 'NoPU_' + version
    elif args.PU > 0:
        version = 'PU{}_'.format(args.PU) + version

    outdir     = args.outdir + '/' + args.process + '_' + version + '/jobs_out'

    time_scale = {'s':1, 'm':60, 'h':60*60, 'd':60*60*24}
    maxRunTime = int(args.maxtime[:-1]) * time_scale[args.maxtime[-1]]
    # mem       = args.memory
    # disk      = args.disk

    mc_frag_dir = args.CMSSW_loc+'/Configuration/GenProduction/python'
    if not os.path.exists(mc_frag_dir):
        os.makedirs(mc_frag_dir)

    mc_frag_name = args.process+'_cfi.py'
    print 'Running process:', args.process
    if not os.path.exists(mc_frag_dir+'/'+mc_frag_name):
        cmd = 'cp '
        cmd += '/eos/user/o/ocerri/BPhysics/MCGeneration/BPH_CMSMCGen/Configuration/GenProduction/python/'
        cmd += mc_frag_name
        cmd += ' '+mc_frag_dir+'/'
        os.system(cmd)
        print 'Compile '+args.CMSSW_loc
        sys.exit()
    else:
        print '--->> I hope you already compiled '+args.CMSSW_loc
        aux = raw_input('Have you? (y/n) ')
        if 'n' in aux:
            exit()

    if not os.path.exists(outdir):
        os.makedirs(outdir)
        os.makedirs(outdir+'/out/')
        os.makedirs(outdir+'/cfg/')
    elif not args.force_production:
        print 'Output dir: "'+outdir+'" exists.'
        aux = raw_input('Continue anyway? (y/n)\n')
        if aux == 'n':
            exit()

    os.system('chmod +x job1023_gen_NoPU_v1.sh')
    print 'Creating submission script\n\n'

    fsub = open('jobs.sub', 'w')
    exec_base = 'executable    = /afs/cern.ch/user/o/ocerri/cernbox/BPhysics/MCGeneration/BPH_CMSMCGen/'
    if args.PU == 0:
        fsub.write(exec_base+'job1023_gen_NoPU_v1.sh')
    elif isinstance(args.PU, int) and args.PU > 0:
        fsub.write(exec_base+'job1023_gen_wPU_v1.sh')
    fsub.write('\n')
    exec_args = str(nev)+' '+str(st_seed)+' $(ProcId) '+args.process+' '+version+' '+args.CMSSW_loc
    if isinstance(args.PU, int) and args.PU > 0:
        exec_args += ' ' + str(args.PU)
    fsub.write('arguments     = ' + exec_args)
    fsub.write('\n')
    fsub.write('output        = {}/out/{}.$(ClusterId).$(ProcId).out'.format(outdir, args.process))
    fsub.write('\n')
    fsub.write('error         = {}/out/{}.$(ClusterId).$(ProcId).err'.format(outdir, args.process))
    fsub.write('\n')
    fsub.write('log           = {}/out/{}.$(ClusterId).$(ProcId).log'.format(outdir, args.process))
    fsub.write('\n')
    fsub.write('+MaxRuntime   = '+str(maxRunTime))
    fsub.write('\n')
    fsub.write('+JobBatchName = '+args.process)
    fsub.write('\n')
    fsub.write('x509userproxy = $ENV(X509_USER_PROXY)')
    fsub.write('\n')
    fsub.write('queue '+str(njobs))
    fsub.write('\n')
    fsub.close()

    output = processCmd('condor_submit jobs.sub')
    os.rename('jobs.sub', outdir+'/cfg/jobs.sub')
