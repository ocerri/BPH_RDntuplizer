import os, sys
from glob import glob
import argparse

parser = argparse.ArgumentParser()

parser.add_argument ('inputDir', type=str, default='tmp/crab_*', help='Input directories', nargs='+')
parser.add_argument ('-t', '--maxTime', type=int, help='Max run time in hours', default=None)
parser.add_argument('--whitelist', type=int, help='Site white list', nargs='+', default=None)
parser.add_argument('--blacklist', type=int, help='Site black list', nargs='+', default=None)


args = parser.parse_args()

for dir in args.inputDir:
    if os.path.isdir(dir):
        print 20*'#' + 50*'-' + 20*'#'
        cmd = 'source /cvmfs/cms.cern.ch/crab3/crab.sh; '
        cmd += 'crab resubmit -d ' + dir
        if not args.maxTime is None:
            cmd += ' --maxjobruntime='+str(60*args.maxTime)
        if not args.whitelist is None:
            cmd += ' --sitewhitelist=' + str(args.whitelist)
        if not args.blacklist is None:
            cmd += ' --siteblacklist=' + str(args.blacklist)
        print cmd
        os.system(cmd)
        print 20*'#' + 50*'-' + 20*'#' + '\n'
