import yaml
import os
import datetime
from glob import glob
import ROOT as rt
from multiprocessing import Pool

def checkFile(floc):
    if not os.path.isfile(floc):
        floc = '/storage/cms' + floc
    if not os.path.isfile(floc):
        # print floc, '(No local copy)'
        return False
    if rt.TFile(floc).IsZombie():
        print floc, '(zombie)'
        return False
    # try:
    #     f = rt.TFile(floc)
    #     tree = f.Get("Events")
    #     for event in tree:
    #         pass
    # except:
    #     print floc, '(could not loop)'

    return True

prod_samples = yaml.load(open('samples.yml'))
for k, d in prod_samples['samples'].iteritems():
    if 'data_' in k:
        continue
        for i in d['parts']:
            if not i in [1,2,4]:
                continue
            dataset = d['dataset'].format(i)
            print '########## {} ##########'.format(dataset)
            fname = 'inputFiles' + dataset.replace('/', '_') + '.txt'
            flist = [l.strip() for l in open(fname, 'r').readlines()]
            pool = Pool(25)
            pool.map(checkFile, flist)
            # nFiles = len(flist)
            # for i_f, f in enumerate(flist):
            #     if i_f % 500 == 0:
            #         print '{}/{}'.format(i_f, nFiles)
            #     checkFile(f)
    else:
        if not k in ['Bd_MuNuDst', 'Bd_DstDd', 'Bu_MuNuDstPi']:
            continue
        dataset = d['dataset']
        print '########## {} ##########'.format(dataset)

        fname = 'inputFiles_' + dataset + '.txt'
        flist = [l.strip() for l in open(fname, 'r').readlines()]
        pool = Pool(25)
        pool.map(checkFile, flist)
        # nFiles = len(flist)
        # for i_f, f in enumerate(flist):
        #     if i_f % 500 == 0:
        #         print '{}/{}'.format(i_f, nFiles)
        #     checkFile(f)
