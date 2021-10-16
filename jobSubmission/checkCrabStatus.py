import os, commands
import argparse
from glob import glob

parser = argparse.ArgumentParser()
parser.add_argument ('inputDir', type=str, default='tmp/crab_*', help='Input dir template for glob', nargs='+')
parser.add_argument ('--report', default=False, action='store_true')
parser.add_argument ('--long', default=False, action='store_true')
parser.add_argument ('--short', default=False, action='store_true')
parser.add_argument ('--verboseErrors', default=False, action='store_true')
parser.add_argument ('--brilcalc', default=False, action='store_true')
args = parser.parse_args()

if args.brilcalc:
    args.report = True

for dir in args.inputDir:
    if os.path.isdir(dir):
        print 20*'#' + 50*'-' + 20*'#'
        cmd = 'source /cvmfs/cms.cern.ch/crab3/crab.sh; '
        cmd += 'crab status -d ' + dir
        if args.verboseErrors:
            cmd += ' --verboseErrors'
        if args.long:
            cmd += ' --long'
        if args.short:
            stats, output = commands.getstatusoutput(cmd)
            outLines = output.split('\n')
            print outLines[0][outLines[0].find('tmp/crab_') + 9: ]
            for i, l in enumerate(outLines):
                if l.startswith('Jobs status:'):
                    jf = 1
                    while '%' in outLines[i+jf]:
                        jf += 1
                    print '\n'.join(outLines[i:i+jf])
                    break
        else:
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
                if os.uname()[1] == 'login-1.hep.caltech.edu':
                    cmd += ' -c web'
                cmd += ' --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json'
                cmd += ' -i ' + dir + '/results/processedLumis.json'
                cmd += ' -o ' + dir + '/results/lumiReport_brilcalc.csv'
                cmd += '; tail ' + dir + '/results/lumiReport_brilcalc.csv'
                os.system(cmd)

                cmd = 'cp ' + dir + '/results/lumiReport_brilcalc.csv'
                if os.uname()[1] == 'login-1.hep.caltech.edu':
                    cmd += ' /storage/user/ocerri/BPhysics/data/cmsRD/lumiReport/'
                else:
                    cmd += ' /afs/cern.ch/user/o/ocerri/cernbox/BPhysics/data/cmsRD/lumiReport/'
                cmd += os.path.basename(dir)[5:] + '_bricalc.csv'
                os.system(cmd)
        if not args.short:
            print 20*'#' + 50*'-' + 20*'#' + '\n\n'
if not args.short:
    print 'For info on exit codes visit: https://twiki.cern.ch/twiki/bin/viewauth/CMSPublic/JobExitCodes'
