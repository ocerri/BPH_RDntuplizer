import os
import argparse
from glob import glob

parser = argparse.ArgumentParser()
parser.add_argument ('inputDir', type=str, default='tmp/crab_*', help='Input dir template for glob')
args = parser.parse_args()

for dir in glob(args.inputDir):
    if os.path.isdir(dir):
        print 20*'#' + 50*'-' + 20*'#'
        cmd = 'source /cvmfs/cms.cern.ch/crab3/crab.sh; '
        cmd += 'crab status -d ' + dir
        cmd += ' --verboseErrors --long'
        os.system(cmd)
        print '\n' + 70*'~' + '\n'
        cmd = 'source /cvmfs/cms.cern.ch/crab3/crab.sh; '
        cmd += 'crab report -d ' + dir
        os.system(cmd)
        print 20*'#' + 50*'-' + 20*'#' + '\n\n'

print 'For info on exit codes visit: https://twiki.cern.ch/twiki/bin/viewauth/CMSPublic/JobExitCodes'
