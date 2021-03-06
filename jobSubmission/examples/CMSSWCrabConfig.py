# To submit: crab submit -c CMSSWCrabConfig.py
from WMCore.Configuration import Configuration
config = Configuration()
config.section_('General')
config.General.transferOutputs = True
config.General.requestName = 'BPH5_Run2018D_RDntuplizer_B02DstMu_190517'
config.section_('JobType')
config.JobType.psetName = '../config/cmssw_cmsRD2018_Tag_B0_MuDmst-pD0bar-kp.py'
config.JobType.pluginName = 'Analysis'
config.JobType.outputFiles = ['B02DstMu_candidates.root']
config.JobType.maxJobRuntimeMin = 240
config.JobType.maxMemoryMB = 2000
config.JobType.allowUndistributedCMSSW = True
config.section_('Data')
#To list the dataset: dasgoclient -query "/ParkingBPH*/*Mar2019*/MINIAOD"
config.Data.inputDataset = '/ParkingBPH5/Run2018D-20Mar2019-v1/MINIAOD'
config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions18/13TeV/ReReco/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt'
config.Data.publication = False
config.Data.unitsPerJob = 5
# config.Data.totalUnits = 3
config.Data.publishDBS = 'https://cmsweb.cern.ch/dbs/prod/phys03/DBSWriter/'
config.Data.splitting = 'FileBased'
config.Data.inputDBS = 'https://cmsweb.cern.ch/dbs/prod/global/DBSReader/'
config.Data.outputDatasetTag = 'BPH5_Run2018D_RDntuplizer_B02DstMu_190517'
config.section_('Site')
# config.Site.blacklist = ['T2_IT_Legnaro']
# config.Site.whitelist = ['T2_IT_Bari']
config.Site.storageSite = 'T2_US_Caltech'
