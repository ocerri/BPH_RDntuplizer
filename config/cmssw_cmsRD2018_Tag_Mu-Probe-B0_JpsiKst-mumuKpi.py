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
    if len(args.inputFiles) == 1:
        with open(args.inputFiles[0]) as f:
            flist = [l for l in f.readlines()]
    else:
        flist = args.inputFiles
else:
    fdefault = os.environ['CMSSW_BASE'] + '/src/ntuplizer/BPH_RDntuplizer/production/'
    fdefault += 'inputFiles_ParkingBPH1_Run2018D-05May2019promptD-v1_MINIAOD.txt'
    with open(fdefault) as f:
        flist = [l[:-1] for l in f.readlines()]
    flist = flist[:10]

print 'Trying to get a local copy'
for i in range(len(flist)):
    if flist[i].startswith('file:'):
        print 'Already set to local'
        continue
    print 'Looking for: /mnt/hadoop' + flist[i]
    if os.path.isfile('/mnt/hadoop' + flist[i]):
        print 'Found'
        flist[i] = 'file:/mnt/hadoop' + flist[i]

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
    outname = 'B2JpsiKst_CAND.root'
elif args.outputFile.startswith('_numEvent'):
    outname = 'B2JpsiKst_CAND' + args.outputFile
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
        verbose = cms.int32(0)
)

process.B2JpsiKstDT = cms.EDProducer("B2JpsiKstDecayTreeProducer",
        trgMuons = cms.InputTag("trgF","trgMuonsMatched", ""),
        verbose = cms.int32(0)
)

process.B2JpsiKstDTFilter = cms.EDFilter("B2JpsiKstDecayTreeFilter",
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
