import os, sys
from glob import glob
import argparse
import humanfriendly

parser = argparse.ArgumentParser()
parser.add_argument ('inputDir', help='list of directories', nargs='+')
parser.add_argument ('-k', '--keepGoing', help='Skip corrupted files', action='store_true')
args = parser.parse_args()
# example: python collectCrabJobs.py /mnt/hadoop/store/user/ocerri/ParkingBPH*/*_RDntuplizer_B2JpsiKst_200123

for dir in args.inputDir:
    print 20*'#' + 50*'-' + 20*'#'
    print 'Direcory:', dir
    if dir[-1] != '/':
        dir += '/'
    dlist = dir.split('/')

    i_BPH = None
    for i,e in enumerate(dlist):
        if e[:-1] == 'ParkingBPH':
            i_BPH = i
    if i_BPH is None:
        print 'No <ParkingBPH> dir found'
        raise

    outpath = '/storage/user/ocerri/BPhysics/data/cmsRD/'
    outpath += dlist[i_BPH]
    if not os.path.isdir(outpath):
        os.makedirs(outpath)

    outname = dlist[i_BPH+1].replace(dlist[i_BPH]+'_', '')
    outpath += '/' + outname

    inputFiles = glob(dir + '*/*/*_CAND_*.root')
    print 'Files to be merged: ', len(inputFiles)
    size_bytes = 0
    for f in inputFiles:
        size_bytes += os.path.getsize(f)
    print 'Total size:', humanfriendly.format_size(size_bytes, binary=True)
    bname = os.path.basename(inputFiles[0])
    lbname = bname.split('_')

    outpath += '_' + lbname[-2] + '.root'
    print 'Output:', outpath

    cmd = 'hadd -f '
    if args.keepGoing:
        cmd += '-k '
    cmd += outpath + ' ' + dir + '*/*/*_CAND_*.root'
    cmd += ' &> merge_' + dlist[i_BPH+1] + '.log &'
    print 'Running:', cmd
    os.system(cmd)
    print 20*'#' + 50*'-' + 20*'#' + '\n\n'
