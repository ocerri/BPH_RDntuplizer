import os
import argparse
from glob import glob

parser = argparse.ArgumentParser()
parser.add_argument ('inputDir', type=str, default='tmp/crab_*', help='Input dir template for glob', nargs='+')
parser.add_argument ('--dry_run', default=False, action='store_true', help='Do not run commands')
parser.add_argument ('--rm', default=False, action='store_true', help='Remove crab directory')
args = parser.parse_args()

for dir in args.inputDir:
    if os.path.isdir(dir):
        print 20*'#' + 50*'-' + 20*'#'
        cmd = 'source /cvmfs/cms.cern.ch/crab3/crab.sh; '
        cmd += 'crab kill -d ' + dir
        if not args.dry_run:
            os.system(cmd)
        else:
            print cmd

        if args.rm and not args.dry_run:
            os.system('rm -rf ' + dir)
        print 20*'#' + 50*'-' + 20*'#' + '\n'
