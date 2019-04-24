import FWCore.ParameterSet.Config as cms
import os

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
#####################   Input    ###################
'''
process.maxEvents = cms.untracked.PSet(
    # input = cms.untracked.int32(50)
    input = cms.untracked.int32(-1)
)

# from glob import glob
# flist = glob('/eos/user/o/ocerri/BPhysics/data/cmsMC_private/BPH_Tag-Bm_D0kpmunu_Probe-B0_MuNuDmst-pD0bar-kp-_NoPU_10-2-3_v1/jobs_out/*MINIAODSIM*.root')
# for i in range(len(flist)):
#     flist[i] = 'file:' + flist[i]
flist =[ 'file:/eos/user/o/ocerri/BPhysics/data/cmsMC_private/BPH_Tag-Bm_D0kpmunu_Probe-B0_MuNuDmst-pD0bar-kp-_NoPU_10-2-3_v1/jobs_out/BPH_Tag-Bm_D0kpmunu_Probe-B0_MuNuDmst-pD0bar-kp-_MINIAODSIM_merged_1-300.root']
process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(tuple(flist)) )

process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')


'''
#####################   Output   ###################
'''

outname = flist[0].replace('jobs_out/', '')
outname = outname.replace('MINIAODSIM', 'BPHRDntuplizer')
# outname = "flat_tree.root"

process.TFileService = cms.Service("TFileService",
      fileName = cms.string(outname),
      closeFileFast = cms.untracked.bool(True)
      )



'''
#################   Sequence    ####################
'''

process.trgBPH = cms.EDProducer("BPHTriggerPathProducer",
        muonCollection = cms.InputTag("slimmedMuons","", "PAT"),
        triggerObjects = cms.InputTag("slimmedPatTrigger"),
        triggerBits = cms.InputTag("TriggerResults","","HLT"),
        muon_charge = cms.int32(-1),
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


process.outA = cms.EDAnalyzer("FlatTreeWriter",
        cmssw = cms.string(cmssw_version),
        verbose = cms.int32(0)

)


process.p = cms.Path(
                    process.trgBPH +
                    process.trgF +
                    process.R2Mmatch +
                    process.R2MmatchFilter +
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
