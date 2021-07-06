import yaml
import os
import datetime
from glob import glob

forceMC = False
prod_samples = yaml.full_load(open('samples.yml'))

for k, d in prod_samples['samples'].iteritems():
    if 'data_' in k:
        for i in d['parts']:
            dataset = d['dataset'].format(i)
            print '########## {} ##########'.format(dataset)
            fname = 'inputFiles' + dataset.replace('/', '_') + '.txt'
            if os.path.isfile(fname):
                print 'Already present'
            else:
                print 'Creating files list'
                cmd = 'dasgoclient --query="file dataset={}" --limit=0 > {}'.format(dataset, fname)
                os.system(cmd)
                print '\nTotal number of files:'
                os.system('cat '+fname+' | wc -l')
            print '\n'
    else:
        dataset = d['dataset']
        print '########## {} ##########'.format(dataset)
        fname = 'inputFiles_' + dataset + '.txt'
        if os.path.isfile(fname):
            if forceMC:
                os.system('rm '+fname)
            else:
                print 'Already present\n\n'
                continue
        for p in d['parts']:
            print p
            if p.endswith('/USER'):
                cmd = 'dasgoclient --query="file dataset={}'.format(p)
                cmd += ' instance=prod/phys03'
                cmd += '"'
                cmd += ' --limit=0 >> ' + fname
                os.system(cmd)
            elif p.endswith('/MINIAODSIM'):
                cmd = 'dasgoclient --query="file dataset={}'.format(p)
                cmd += '"'
                cmd += ' --limit=0 >> ' + fname
                os.system(cmd)
            else:
                with open(fname, 'a') as f:
                    for name in glob(p):
                        f.write(name+'\n')

        print '\nTotal number of files:'
        os.system('cat '+fname+' | wc -l')
        print '\n'
