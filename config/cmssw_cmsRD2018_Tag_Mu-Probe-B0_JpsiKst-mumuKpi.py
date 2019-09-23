import os, sys
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

from Configuration.StandardSequences.Eras import eras
process = cms.Process('BPHRDntuplizer', eras.Run2_2018)
cmssw_version = os.environ['CMSSW_VERSION']
# import of standard configurations
process.load('FWCore.MessageService.MessageLogger_cfi')


process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '102X_dataRun2_v11', '')

'''
############ Command line args ################
'''

args = VarParsing.VarParsing('analysis')
args.register('inputFile', '', args.multiplicity.list, args.varType.string, "Input file or template for glob")
args.outputFile = ''
args.parseArguments()


'''
#####################   Input    ###################
'''
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(args.maxEvents)
)

from glob import glob
if args.inputFile:
    flist = args.inputFile
elif args.inputFiles:
    if len(args.inputFiles) == 1:
        with open(args.inputFiles[0]) as f:
            flist = [l for l in f.readlines()]
    else:
        flist = args.inputFiles
else:
    flist = glob('/eos/cms/store/data/Run2018C/ParkingBPH1/MINIAOD/05May2019-v1/*/*.root')
    for i in range(len(flist)):
        flist[i] = 'file:' + flist[i]

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(tuple(flist))
                           )

# import FWCore.PythonUtilities.LumiList as LumiList
# process.source.lumisToProcess = LumiList.LumiList(filename = '/afs/cern.ch/user/o/ocerri/work/CMSSW_10_2_3/src/ntuplizer/BPH_RDntuplizer/production/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt').getVLuminosityBlockRange()


'''
#####################   Output   ###################
'''
if args.outputFile == '.root':
    outname = 'B2JpsiKst_CAND.root'
else:
    outname = args.outputFile

process.TFileService = cms.Service("TFileService",
      fileName = cms.string(outname),
      closeFileFast = cms.untracked.bool(True)
      )



'''
#################   Sequence    ####################
'''

process.trgBPH = cms.EDProducer("BPHTriggerPathProducer",
        muonCollection = cms.InputTag("slimmedMuons","", ""),
        vertexCollection = cms.InputTag("offlineSlimmedPrimaryVertices","", ""),
        triggerObjects = cms.InputTag("slimmedPatTrigger"),
        triggerBits = cms.InputTag("TriggerResults","","HLT"),
        muon_charge = cms.int32(0),
        verbose = cms.int32(0)
)

process.trgF = cms.EDFilter("BPHTriggerPathFilter",
        trgMuons = cms.InputTag("trgBPH","trgMuonsMatched", "")
)

process.B2JpsiKstDT = cms.EDProducer("B2JpsiKstDecayTreeProducer",
        trgMuons = cms.InputTag("trgBPH","trgMuonsMatched", ""),
        verbose = cms.int32(1)
)

process.B2JpsiKstDTFilter = cms.EDFilter("B2JpsiKstDecayTreeFilter",
        verbose = cms.int32(1)
)

cfg_name = os.path.basename(sys.argv[0])
f = open(os.environ['CMSSW_BASE']+'/src/ntuplizer/BPH_RDntuplizer/.git/logs/HEAD')
commit_hash = f.readlines()[-1][:-1].split(' ')[1]
process.outA = cms.EDAnalyzer("FlatTreeWriter",
        cmssw = cms.string(os.environ['CMSSW_VERSION']),
        cfg_name = cms.string(cfg_name),
        commit_hash = cms.string(commit_hash),
        verbose = cms.int32(1)
)


process.p = cms.Path(
                    process.trgBPH +
                    process.trgF +
                    process.B2JpsiKstDT +
                    process.B2JpsiKstDTFilter +
                    process.outA
                    )

# DEBUG -- dump the event content
# process.output = cms.OutputModule(
#                 "PoolOutputModule",
#                       fileName = cms.untracked.string('edm_output.root'),
#                       )
# process.output_step = cms.EndPath(process.output)
#
# process.schedule = cms.Schedule(
# 		process.p,
# 		process.output_step)


'''
#############   Overall settings    ################
'''

process.MessageLogger.cerr.FwkReport.reportEvery = 1000
