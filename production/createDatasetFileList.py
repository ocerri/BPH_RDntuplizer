import yaml
import os
import datetime

prod_samples = yaml.load(open('samples.yml'))

for k, d in prod_samples['samples'].iteritems():
    if 'data_' in k:
        for i in d['parts']:
            dataset = d['dataset'].format(i)
            print '########## {} ##########\n'.format(dataset)
            fname = 'inputFiles' + dataset.replace('/', '_')
            cmd = 'das_client --query="file dataset={}" --limit=0 > {}.txt'.format(dataset, fname)
            os.system(cmd)
            print '\n'
