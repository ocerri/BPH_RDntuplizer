import os
import argparse
from glob import glob

parser = argparse.ArgumentParser()
parser.add_argument ('inputDir', type=str, default='tmp/crab_*', help='Input dir template for glob')
parser.add_argument ('--report', default=False, action='store_true')
parser.add_argument ('--long', default=False, action='store_true')
parser.add_argument ('--brilcalc', default=False, action='store_true')
args = parser.parse_args()

if args.brilcalc:
    args.report = True

for dir in glob(args.inputDir):
    if os.path.isdir(dir):
        print 20*'#' + 50*'-' + 20*'#'
        cmd = 'source /cvmfs/cms.cern.ch/crab3/crab.sh; '
        cmd += 'crab status -d ' + dir
        cmd += ' --verboseErrors'
        if args.long:
            cmd += ' --long'
        os.system(cmd)

        if args.report:
            print '\n' + 70*'~' + '\n'
            cmd = 'source /cvmfs/cms.cern.ch/crab3/crab.sh; '
            cmd += 'crab report -d ' + dir
            os.system(cmd)
            if args.brilcalc:
                print 'Running brilcalc...\n'
                cmd = 'export PATH=$HOME/.local/bin:/cvmfs/cms-bril.cern.ch/brilconda/bin:$PATH'
                cmd += '; brilcalc lumi -u /fb'
                cmd += ' --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json'
                cmd += ' -i ' + dir + '/results/processedLumis.json'
                cmd += ' -o ' + dir + '/results/lumireport_brilcalc.csv'
                cmd += '; tail ' + dir + '/results/lumireport_brilcalc.csv'
                os.system(cmd)
        print 20*'#' + 50*'-' + 20*'#' + '\n\n'

print 'For info on exit codes visit: https://twiki.cern.ch/twiki/bin/viewauth/CMSPublic/JobExitCodes'
