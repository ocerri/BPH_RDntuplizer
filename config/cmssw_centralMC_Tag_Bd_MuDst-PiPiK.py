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

    # fdefault += 'inputFiles_CP_BdToDstarMuNu_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen.txt'
    # fdefault += 'inputFiles_CP_BdToDstarTauNu_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen.txt'
    # fdefault += 'inputFiles_CP_BdToDstDs_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen.txt'
    # fdefault += 'inputFiles_CP_BuToMuNuDstPi_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen.txt'
    fdefault += 'inputFiles_CP_BdToMuNuDstPiPi_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen_v3.txt'

    with open(fdefault) as f:
        flist = [l[:-1] for l in f.readlines()]
    # flist = flist[:5]
    flist = flist[-5:]

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
    outname = 'B2DstMu_CAND.root'
elif args.outputFile.startswith('_numEvent'):
    outname = 'B2DstMu_CAND' + args.outputFile
else:
    outname = args.outputFile

process.TFileService = cms.Service("TFileService",
      fileName = cms.string(outname),
      closeFileFast = cms.untracked.bool(True)
      )



'''
#################   Sequence    ####################
'''

process.trgF = cms.EDFilter("TriggerMuonsFilter",
        muon_charge = cms.int32(0),
        isMC = cms.int32(1),
        verbose = cms.int32(0)
)

process.B2MuDstDT = cms.EDProducer("B2DstMuDecayTreeProducer",
        trgMuons = cms.InputTag("trgF","trgMuonsMatched", ""),
        charge_muon = cms.int32(0), # if 0 accept both charges
        charge_K = cms.int32(+1), # charges relative to the muon charge
        charge_pi = cms.int32(-1), # charges relative to the muon charge
        charge_pis = cms.int32(-1), # charges relative to the muon charge
        verbose = cms.int32(0)
)

process.B2MuDstDTFilter = cms.EDFilter("B2DstMuDecayTreeFilter",
        verbose = cms.int32(0)
)

process.MCpart = cms.EDProducer("MCTruthB2DstMuProducer",
        decayTreeVecOut = cms.InputTag("B2MuDstDT","outputVecNtuplizer", ""),
        verbose = cms.int32(0)
)

process.HammerWeights = cms.EDProducer("HammerWeightsProducer",
        decayOfInterest = cms.vstring('BD*MuNu', 'BD*TauNu', 'BD**MuNu', 'BD2sallMuNu'),
        inputFFScheme = cms.vstring(#'BD', 'ISGW2',
                                    'BD*', 'ISGW2',
                                    # 'BD**', 'ISGW2',
                                    'BD**1', 'ISGW2',
                                    'BD**2*', 'ISGW2',
                                    'BD**0*', 'ISGW2',
                                    'BD**1*', 'ISGW2',
                                    'BD(2s)', 'ISGW2',
                                    'BD(2s)*', 'ISGW2',
        ),
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
                    process.trgF +
                    process.B2MuDstDT +
                    process.B2MuDstDTFilter+
                    process.MCpart +
                    process.HammerWeights +
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
# process.MessageLogger.cerr.threshold = "DEBUG"
# process.MessageLogger.debugModules = ["*"]
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
