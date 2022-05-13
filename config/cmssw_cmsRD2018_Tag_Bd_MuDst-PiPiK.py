import os, sys
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

from Configuration.StandardSequences.Eras import eras
process = cms.Process('BPHRDntuplizer', eras.Run2_2018)
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
args.register('useLocalLumiList', 1, args.multiplicity.singleton, args.varType.int, "Flag to use local lumi list")
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
    if len(args.inputFiles) == 1 and args.inputFiles[0].endswith('.txt'):
        with open(args.inputFiles[0]) as f:
            flist = [l[:-1] for l in f.readlines()]
    else:
        flist = args.inputFiles
else:
    fdefault = os.environ['CMSSW_BASE'] + '/src/ntuplizer/BPH_RDntuplizer/production/'
    fdefault += 'inputFiles_ParkingBPH1_Run2018D-05May2019promptD-v1_MINIAOD.txt'
    with open(fdefault) as f:
        flist = [l[:-1] for l in f.readlines()]
    flist = flist[:2]
    # flist = ['/store/data/Run2018D/ParkingBPH5/MINIAOD/05May2019promptD-v1/230005/F31C64F8-800A-974E-8458-FC06FEC9442D.root']

# print 'Trying to get a local copy'
for i in range(len(flist)):
    if os.path.isfile(flist[i]):
        flist[i] = 'file:' + flist[i]
    elif flist[i].startswith('/store/data/'):
        if os.path.isfile('/storage/cms' + flist[i]):
            flist[i] = 'file:' + '/storage/cms' + flist[i]

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(tuple(flist))
                           )

if args.useLocalLumiList:
    lumiListFile = os.environ['CMSSW_BASE'] + '/src/ntuplizer/BPH_RDntuplizer/production/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt'
    if os.path.isfile(lumiListFile):
        import FWCore.PythonUtilities.LumiList as LumiList
        process.source.lumisToProcess = LumiList.LumiList(filename = lumiListFile).getVLuminosityBlockRange()

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
        isMC = cms.int32(0),
        verbose = cms.int32(0)
)


process.B2MuDstDT = cms.EDProducer("B2DstMuDecayTreeProducer",
        trgMuons = cms.InputTag("trgF","trgMuonsMatched", ""),
        charge_muon = cms.int32(0), # if 0 accept both charges
        charge_K = cms.int32(+1), # charges relative to the muon charge
        charge_pi = cms.int32(-1), # charges relative to the muon charge
        charge_pis = cms.int32(-1), # charges relative to the muon charge
        verbose = cms.int32(0),
        isMC = cms.int32(0),
)

process.B2MuDstDTFilter = cms.EDFilter("B2DstMuDecayTreeFilter",
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

process.Timing = cms.Service("Timing",
  summaryOnly = cms.untracked.bool(True),
  useJobReport = cms.untracked.bool(True)
)

process.options = cms.untracked.PSet(
 wantSummary = cms.untracked.bool(True),
 numberOfThreads = cms.untracked.uint32(1),
  numberOfStreams = cms.untracked.uint32(1),
)
