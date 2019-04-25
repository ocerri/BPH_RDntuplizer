import os
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

from Configuration.StandardSequences.Eras import eras
process = cms.Process('BPHRDntuplizer', eras.Run2_2018)
cmssw_version = os.environ['CMSSW_VERSION']
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
    # input = cms.untracked.int32(50)
    input = cms.untracked.int32(-1)
)

from glob import glob
if args.inputFile:
    flist = args.inputFile
if args.inputFiles:
    if len(args.inputFiles) == 1:
        with open(args.inputFiles[0]) as f:
            flist = [l for l in f.readlines()]
    else:
        flist = args.inputFiles
else:
    flist = glob('/eos/user/o/ocerri/BPhysics/data/cmsMC_private/BPH_Tag-B0_MuNuDmst-pD0bar-kp_13TeV-pythia8_SoftQCD_PTFilter5_0p0-evtgen_HQET2_central_PU35_10-2-3_v0/jobs_out/*MINIAODSIM*.root')

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
    outname = '/eos/user/o/ocerri/BPhysics/data/cmsMC_private/BPH_Tag-B0_MuNuDmst-pD0bar-kp_13TeV-pythia8_SoftQCD_PTFilter5_0p0-evtgen_HQET2_central_PU35_10-2-3_v0/B02DstMu_candidates.root'
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
        muon_charge = cms.int32(1),
        verbose = cms.int32(0)
)

process.trgF = cms.EDFilter("BPHTriggerPathFilter",
        trgMuons = cms.InputTag("trgBPH","trgMuonsMatched", "")
)


process.R2Mmatch = cms.EDProducer("RECOMCmatchDecayRecoProducer",
        verbose = cms.int32(0)
)

process.R2MmatchFilter = cms.EDFilter("RECOMCmatchDecayRecoFilter",
        verbose = cms.int32(0)
)

process.B2MuDstDT = cms.EDProducer("B2DstMuDecayTreeProducer",
        trgMuons = cms.InputTag("trgBPH","trgMuonsMatched", ""),
        verbose = cms.int32(0)
)

process.B2MuDstDTFilter = cms.EDFilter("B2DstMuDecayTreeFilter", verbose = cms.int32(1))

process.outA = cms.EDAnalyzer("FlatTreeWriter",
        cmssw = cms.string(cmssw_version),
        verbose = cms.int32(0)

)


process.p = cms.Path(
                    process.trgBPH +
                    process.trgF +
                    # process.R2Mmatch +
                    # process.R2MmatchFilter +
                    process.B2MuDstDT +
                    process.B2MuDstDTFilter+
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
