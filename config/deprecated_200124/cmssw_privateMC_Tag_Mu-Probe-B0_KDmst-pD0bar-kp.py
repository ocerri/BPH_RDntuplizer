import os, sys
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

from Configuration.StandardSequences.Eras import eras
process = cms.Process('BPHRDntuplizer', eras.Run2_2018)
# import of standard configurations
process.load('FWCore.MessageService.MessageLogger_cfi')


# Needed for transient track builder
# process.load('Configuration.StandardSequences.Services_cff')
# process.load('Configuration.EventContent.EventContent_cff')
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
# process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '102X_upgrade2018_realistic_v12', '')

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
    if len(args.inputFile) == 1 and '*' in args.inputFile[0]:
        flist = glob(args.inputFile[0])
    else:
        flist = args.inputFile
elif args.inputFiles:
    if len(args.inputFiles) == 1:
        with open(args.inputFiles[0]) as f:
            flist = [l for l in f.readlines()]
    else:
        flist = args.inputFiles
else:
    flist = glob('/afs/cern.ch/user/o/ocerri/cernbox/BPhysics/data/cmsMC_private/BPH_Tag-Mu_Probe-B0_KDmst-pD0bar-kp_13TeV-pythia8_Hardbbbar_PTFilter5_0p0-evtgen_SVS_PU0_10-2-3/*MINIAODSIM*.root')

for i in range(len(flist)):
    flist[i] = 'file:' + flist[i]

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(tuple(flist)),
                            inputCommands=cms.untracked.vstring('keep *',
                                                                'drop GenLumiInfoHeader_generator__SIM')
                           )
process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')


'''
#####################   Output   ###################
'''

if args.outputFile == '.root':
    outname = 'B2DstK_CAND.root'
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

process.B2DstKDT = cms.EDProducer("B2DstKDecayTreeProducer",
        trgMuons = cms.InputTag("trgBPH","trgMuonsMatched", ""),
        verbose = cms.int32(0)
)

process.B2DstKDTFilter = cms.EDFilter("B2DstKDecayTreeFilter",
        verbose = cms.int32(0)
)

process.MCpart = cms.EDProducer("MCTruthB2DstKProducer",
        trgMuons = cms.InputTag("trgBPH","trgMuonsMatched", ""),
        verbose = cms.int32(0)
)

cfg_name = os.path.basename(sys.argv[0])
f = open(os.environ['CMSSW_BASE']+'/src/ntuplizer/BPH_RDntuplizer/.git/logs/HEAD')
commit_hash = f.readlines()[-1][:-1].split(' ')[1]
process.outA = cms.EDAnalyzer("FlatTreeWriter",
        cmssw = cms.string(os.environ['CMSSW_VERSION']),
        cfg_name = cms.string(cfg_name),
        commit_hash = cms.string(commit_hash),
        verbose = cms.int32(0)
)


process.p = cms.Path(
                    process.trgBPH +
                    process.trgF +
                    process.B2DstKDT +
                    process.B2DstKDTFilter +
                    process.MCpart +
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
