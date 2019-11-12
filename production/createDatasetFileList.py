import yaml
import os
import datetime

prod_samples = yaml.load(open('samples.yml'))

for k, d in prod_samples['samples'].iteritems():
    if 'data_' in k:
        continue
        for i in d['parts']:
            dataset = d['dataset'].format(i)
            print '########## {} ##########\n'.format(dataset)
            fname = 'inputFiles' + dataset.replace('/', '_') + '.txt'
            cmd = 'das_client --query="file dataset={}" --limit=0 > {}'.format(dataset, fname)
            os.system(cmd)
    else:
        dataset = d['dataset']
        print '########## {} ##########\n'.format(dataset)
        fname = 'inputFiles_' + dataset + '.txt'
        if os.path.isfile(fname):
            os.system('rm '+fname)
        for p in d['parts']:
            print p
            if p.endswith('/USER'):
                cmd = 'das_client --query="file dataset={}'.format(p)
                cmd += ' instance=prod/phys03'
                cmd += '"'
                cmd += ' --limit=0 >> ' + fname
                os.system(cmd)
            else:
                os.system('ls ' + p + ' >> ' + fname)

    print '\nTotal number of files:'
    os.system('cat '+fname+' | wc -l')
    print '\n'
