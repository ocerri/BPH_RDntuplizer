import yaml
import os
import datetime

# tag = 'B2DstMu'
# tag = 'B2DstK'
tag = 'B2JpsiKst'

cfg = {'B2DstMu': 'cmssw_cmsRD2018_Tag_B0_MuDmst-pD0bar-kp.py',
       'B2DstK': 'cmssw_cmsRD2018_Tag_Mu-Probe-B0_KDmst-pD0bar-kp.py',
       'B2JpsiKst': 'cmssw_cmsRD2018_Tag_Mu-Probe-B0_JpsiKst-mumuKpi.py'
       }

date = datetime.datetime.today()
date_str = '{}{:02}{:02}'.format(date.year%100, date.month, date.day)

prod_samples = yaml.load(open('../production/samples.yml'))

def dumpCrabConfig(k, dataset, settings):
    ds_list = dataset.split('/')
    fname = 'cfg' + dataset.replace('/', '_') + '.py'
    fout = open('tmp/' + fname, 'w')
    fout.write(
'''# To submit: crab submit -c CMSSWCrabConfig.py
from WMCore.Configuration import Configuration
config = Configuration()
config.section_('General')
config.General.transferOutputs = True
config.General.transferLogs = True
'''
    )
    dataset_tag = '{}_{}_RDntuplizer_{}_{}'.format(ds_list[1], ds_list[2], tag, date_str)
    fout.write("config.General.requestName = '{}'".format(dataset_tag))
    fout.write('\n')
    fout.write("config.section_('JobType')")
    fout.write('\n')
    fout.write("config.JobType.psetName = '{}/src/ntuplizer/BPH_RDntuplizer/config/{}'".format(os.environ['CMSSW_BASE'], cfg[tag]))
    fout.write('\n')
    fout.write("config.JobType.pluginName = 'Analysis'")
    fout.write('\n')
    fout.write("config.JobType.outputFiles = ['{}_CAND.root']".format(tag))
    # fout.write('\n')
    # fout.write("config.JobType.maxJobRuntimeMin = {}".format(60*int(settings['common']['data']['maxRunTime'])))
    fout.write('\n')
    fout.write("config.JobType.maxMemoryMB = 2000")
    fout.write('\n')
    fout.write("config.JobType.allowUndistributedCMSSW = True")
    fout.write('\n')
    fout.write("config.section_('Data')")
    fout.write('\n')
    fout.write("config.Data.inputDataset = '{}'".format(dataset))
    fout.write('\n')
    fout.write("config.Data.lumiMask = '{}'".format(settings['common']['data']['lumimask']))
    fout.write('\n')
    fout.write("config.Data.publication = False")
    fout.write('\n')
    fout.write("config.Data.unitsPerJob = {}".format(60*int(settings['common']['data']['unitsPerJob'])))
    fout.write('\n')
    fout.write("config.Data.publishDBS = 'https://cmsweb.cern.ch/dbs/prod/phys03/DBSWriter/'")
    fout.write('\n')
    fout.write("config.Data.splitting = 'Automatic'")
    fout.write('\n')
    fout.write("config.Data.inputDBS = 'https://cmsweb.cern.ch/dbs/prod/global/DBSReader/'")
    fout.write('\n')
    fout.write("config.Data.outputDatasetTag = '{}'".format(dataset_tag))
    fout.write('\n')
    fout.write("config.Data.ignoreLocality = True") # requires a whitelist
    fout.write('\n')
    fout.write("config.Data.allowNonValidInputDataset = True")
    fout.write('\n')
    fout.write("config.section_('Site')")
    fout.write('\n')
    fout.write("config.Site.storageSite = 'T2_US_Caltech'")
    fout.write('\n')
    fout.write("config.Site.whitelist = ['T2_US_*']")
    fout.write('\n')
    fout.close()

    return fname, dataset_tag

if not os.path.isdir('tmp'):
    os.system('mkdir tmp')

for k, d in prod_samples['samples'].iteritems():
    if 'data_Run2018' in k:
        for i in d['parts']:
            if int(i) != 1: continue
            dataset = d['dataset'].format(i)
            print '\n########## {} ##########\n'.format(dataset)

            fname, dset_tag = dumpCrabConfig(k, dataset, prod_samples)
            if os.path.isdir('tmp/crab_' + dset_tag):
                os.system('rm -rf tmp/crab_' + dset_tag)
            print '---> Config dump'
            os.system('cat tmp/' + fname)
            print ''
            cmd = 'cd tmp; source /cvmfs/cms.cern.ch/crab3/crab.sh;'
            cmd += ' crab submit -c ' + fname
            os.system(cmd)
            print '\n\n'
