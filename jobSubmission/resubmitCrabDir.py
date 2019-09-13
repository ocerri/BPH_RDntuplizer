import os, sys
from glob import glob
import argparse

parser = argparse.ArgumentParser()

parser.add_argument ('-i', '--inputDir', type=str, default='tmp/crab_*', help='Input dir template for glob')
parser.add_argument ('-t', '--maxTime', type=int, help='Max run time in hours', default=None)

args = parser.parse_args()

for dir in glob(args.inputDir):
    if os.path.isdir(dir):
        print 20*'#' + 50*'-' + 20*'#'
        cmd = 'source /cvmfs/cms.cern.ch/crab3/crab.sh; '
        cmd += 'crab resubmit -d ' + dir
        if not args.maxTime is None:
            cmd += '--maxjobruntime='+str(60*args.maxTime)
        print cmd
        os.system(cmd)
        print 20*'#' + 50*'-' + 20*'#' + '\n'
