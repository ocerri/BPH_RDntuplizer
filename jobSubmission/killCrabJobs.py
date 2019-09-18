import os
import argparse
from glob import glob

parser = argparse.ArgumentParser()
parser.add_argument ('inputDir', type=str, default='tmp/crab_*', help='Input dir template for glob')
parser.add_argument ('--rm', default=False, action='store_true', help='Remove crab directoru')
args = parser.parse_args()

for dir in glob(args.inputDir):
    if os.path.isdir(dir):
        print 20*'#' + 50*'-' + 20*'#'
        cmd = 'source /cvmfs/cms.cern.ch/crab3/crab.sh; '
        cmd += 'crab kill -d ' + dir
        os.system(cmd)

        if args.rm:
            os.system('rm -rf ' + dir)
        print 20*'#' + 50*'-' + 20*'#' + '\n'

print 'For info on exit codes visit: https://twiki.cern.ch/twiki/bin/viewauth/CMSPublic/JobExitCodes'
