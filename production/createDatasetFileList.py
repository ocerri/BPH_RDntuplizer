import yaml
import os
import datetime
from glob import glob

forceMC = False
prod_samples = yaml.load(open('samples.yml'))

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
                cmd = 'das_client --query="file dataset={}" --limit=0 > {}'.format(dataset, fname)
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
                cmd = 'das_client --query="file dataset={}'.format(p)
                cmd += ' instance=prod/phys03'
                cmd += '"'
                cmd += ' --limit=0 >> ' + fname
                os.system(cmd)
            else:
                with open(fname, 'a') as f:
                    for name in glob(p):
                        f.write(name+'\n')

            if p == '/cmsMC_private/ocerri-BPH_Tag-Probe_B0_JpsiKst-mumuKpi-kp_13TeV-pythia8_Hardbbbar_PTFilter5_0p0-evtgen_SVV-c21dec93027231dc6f615dfe5c662834/USER':
                print 'Removing by hand the 191005 part of the JPsiKst dataset'
                with open(fname, 'r') as f:
                    lines = f.readlines()
                with open(fname, 'w') as f:
                    for line in lines:
                        if not '/191005_' in line:
                            f.write(line)
            elif 'ocerri-BPH_Tag-B0_TauNuDmst' in p:
                print "Removing corrupted files"
                removed = 0
                with open(fname, 'r') as f:
                    lines = f.readlines()
                with open(fname, 'w') as f:
                    for line in lines:
                        condition = not '/191111_232736/0001/out_MINIAODSIM_1515.root' in line
                        condition *= not '191111_232736/0003/out_MINIAODSIM_3049.root' in line
                        condition *= not '191205_195222/0005/out_MINIAODSIM_5429.root' in line
                        condition *= not '191205_195222/0002/out_MINIAODSIM_2764.root' in line
                        if condition:
                            f.write(line)
                        else: removed += 1
                print 'Removed', removed
            elif 'ocerri-BPH_Tag-B0_MuNuDmst-pD0bar-kp_13TeV-pythia8_Hardbbbar_PTFilter5_0p0-evtgen_HQET2_central' in p:
                print "Removing corrupted files"
                with open(fname, 'r') as f:
                    lines = f.readlines()
                with open(fname, 'w') as f:
                    for line in lines:
                        condition = not '191119_023455/0001/out_MINIAODSIM_1273.root' in line
                        condition *= not '191119_023455/0001/out_MINIAODSIM_1082.root' in line
                        condition *= not '200225_164158/0002/out_MINIAODSIM_2278.root' in line
                        condition *= not '200225_164158/0008/out_MINIAODSIM_8661.root' in line
                        condition *= not '200225_164158/0008/out_MINIAODSIM_8664.root' in line
                        if condition:
                            f.write(line)
            elif 'ocerri-BPH_Tag-B0_DmstHc-pD0bar-kp-Hc2mu_13TeV-pythia8_Hardbbbar' in p:
                print "Removing corrupted files"
                with open(fname, 'r') as f:
                    lines = f.readlines()
                with open(fname, 'w') as f:
                    for line in lines:
                        condition = not '200118_025132/0002/out_MINIAODSIM_2393.root' in line
                        if condition:
                            f.write(line)

        print '\nTotal number of files:'
        os.system('cat '+fname+' | wc -l')
        print '\n'
