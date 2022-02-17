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
    if len(args.inputFiles) == 1 and args.inputFiles[0].endswith('.txt'):
        with open(args.inputFiles[0]) as f:
            flist = [l[:-1] for l in f.readlines()]
    else:
        flist = args.inputFiles
else:
    fdefault = os.environ['CMSSW_BASE'] + '/src/ntuplizer/BPH_RDntuplizer/production/'
    # fdefault += 'inputFiles_BP_Tag_B0_MuNuDmst_Hardbbbar_evtgen_ISGW2_PUc0_10-2-3.txt'
    fdefault += 'inputFiles_CP_General_BdToJpsiKstar_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen.txt'
    with open(fdefault) as f:
        flist = [l[:-1] for l in f.readlines()]
    flist = flist[:3]

for i in range(len(flist)):
    if os.path.isfile(flist[i]):
        flist[i] = 'file:' + flist[i]
    elif flist[i].startswith('/store/mc/'):
        if os.path.isfile('/storage/cms' + flist[i]):
            flist[i] = 'file:' + '/storage/cms' + flist[i]

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(tuple(flist)),
                            # inputCommands=cms.untracked.vstring('keep *',
                            #                                     'drop GenLumiInfoHeader_generator__SIM'),
                            # skipBadFiles=cms.untracked.bool(True)
                           )
# process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')


'''
#####################   Output   ###################
'''
if args.outputFile == '.root':
    outname = 'TagAndProbeTrigger_CAND.root'
elif args.outputFile.startswith('_numEvent'):
    outname = 'TagAndProbeTrigger_CAND' + args.outputFile
else:
    outname = args.outputFile

process.TFileService = cms.Service("TFileService",
      fileName = cms.string(outname),
      closeFileFast = cms.untracked.bool(True)
      )



'''
#################   Sequence    ####################
'''

process.l1bits=cms.EDProducer("L1TriggerResultsConverter",
                              src=cms.InputTag("gtStage2Digis"),
                              # src_ext=cms.InputTag("gtStage2Digis"),
                              # storeUnprefireableBit=cms.bool(True),
                              legacyL1=cms.bool(False),
)

process.TnP = cms.EDFilter("TagAndProbeTriggerMuonFilter",
        muonIDScaleFactors = cms.int32(1),
        requireTag = cms.int32(1),
        isMC = cms.int32(1),
        verbose = cms.int32(0)
)


process.p = cms.Path(
                    process.l1bits +
                    process.TnP
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

process.MessageLogger.cerr.FwkReport.reportEvery = 10000
