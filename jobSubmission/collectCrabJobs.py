import os, sys
from glob import glob
import argparse
import humanfriendly

parser = argparse.ArgumentParser()

parser.add_argument ('inputDir', help='Input dir template for glob', nargs='+')

args = parser.parse_args()

for dir in args.inputDir:
    print 20*'#' + 50*'-' + 20*'#'
    print 'Direcory:', dir
    dlist = dir.split('/')

    i_BPH = None
    for i,e in enumerate(dlist):
        if e[:-1] == 'ParkingBPH':
            i_BPH = i
    if i_BPH is None:
        print 'No <ParkingBPH> dir found'
        raise

    outpath = '/storage/user/ocerri/cmsRD/'
    outpath += dlist[i_BPH]
    if not os.path.isdir(outpath):
        os.makedirs(outpath)

    outname = dlist[i_BPH+1].replace(dlist[i_BPH]+'_', '')
    outpath += '/' + outname

    inputFiles = glob(dir + '*/*/*.root')
    print 'Files to be merged: ', len(inputFiles)
    size_bytes = 0
    for f in inputFiles:
        size_bytes += os.path.getsize(f)
    print 'Total size:', humanfriendly.format_size(size_bytes, binary=True)
    bname = os.path.basename(inputFiles[0])
    lbname = bname.split('_')

    outpath += '_' + lbname[-2] + '.root'
    print 'Output:', outpath

    cmd = 'hadd ' + outpath + ' ' + dir + '*/*/*.root'
    cmd += ' &> merge_' + dlist[i_BPH+1] + '.log &'
    print 'Running:', cmd
    os.system(cmd)
    print 20*'#' + 50*'-' + 20*'#' + '\n\n'
