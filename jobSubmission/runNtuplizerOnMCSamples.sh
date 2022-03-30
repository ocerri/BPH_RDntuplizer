#!/usr/bin/env bash
# Script to submit jobs to condor for running the ntuplizer. To run it just
# run:
#
#     # sh ./jobSubmission/runNtuplizerOnMCSamples.sh
#
# Make sure you have updated the input files by running:
#
#     # cd production
#     # python createDatasetFileList.py
#
# first.

outLoc=/storage/af/group/rdst_analysis/BPhysics/data/cmsMC

# Trigger tag and probe ntuples
# config=config/cmssw_centralMC_TagAndProbeTrigger.py
# ntuplesName=ntuples_TagAndProbeTrigger_220217

# Bd -> Jpsi K* ntuples
# config=config/cmssw_centralMC_Bd_JpsiKst-mumuKpi.py
# ntuplesName=ntuples_Bd2JpsiKst_220328

# Main R(D*) analysis ntuples
config=config/cmssw_centralMC_Tag_Bd_MuDst-PiPiK.py
ntuplesName=ntuples_B2DstMu_220326_wBGL_v0

# nFilesPerJob=100
# maxTime=120m

nFilesPerJob=8
maxTime=48h
declare -a processes=(
    # Ancillary measurments samples --> Should be run N = 3
    # "CP_General_BdToJpsiKstar_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen"
    # "CP_General_BuToJpsiK_BMuonFilter_TuneCP5_13TeV-pythia8-evtgen"
    # "CP_General_MuEnriched_HardQCDall_TuneCP5_13TeV-pythia8"
    #
    # Central production --> Should be run N = 3
    "CP_BdToDstarMuNu_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen"
    "CP_BdToDstarTauNu_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen"
    # "CP_BuToMuNuDstPi_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen"
    # "CP_BdToMuNuDstPi_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen"
    # "CP_BuToTauNuDstPi_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen"
    # "CP_BdToTauNuDstPi_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen"
    # "CP_BsToMuNuDstK_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen"
    # "CP_BsToTauNuDstK_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen"
    # "CP_BdToDstDu_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen"
    # "CP_BdToDstDd_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen"
    # "CP_BdToDstDs_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen"
    # "CP_BuToDstDu_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen"
    # "CP_BuToDstDd_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen"
    # "CP_BsToDstDs_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen"
    #
    # Private production --> Should be run N = 100
    # "CP_BdToMuNuDstPiPi_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen_v3"
    # "CP_BuToMuNuDstPiPi_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen_v3"
    # "CP_BdToTauNuDstPiPi_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen"
    # "CP_BuToTauNuDstPiPi_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen"
    # "BParking_Tag_Bd_DDs1_SoftQCDnonD_TuneCP5_13TeV-pythia8"
    # "BParking_Tag_Bu_DDs1_SoftQCDnonD_TuneCP5_13TeV-pythia8"
    # "BParking_Tag_B_DstDXX_SoftQCDnonD_TuneCP5_13TeV-pythia8"
)

for process in "${processes[@]}"; do
    echo $process
    output_dir=$outLoc/$process/$ntuplesName
    mkdir -p $output_dir
    python jobSubmission/create-condor-jobs -i production/inputFiles_$process.txt -o $output_dir/out_CAND.root -c $config -t $ntuplesName -N $nFilesPerJob --maxtime $maxTime
    sleep 1
done


# To actually submit the jobs run: ./jobSubmission/submit-condor-jobs --max-jobs 1000
# To check the status of the submission:
#   sqlite3 ~/state.db
#   sqlite3 -column -header ~/state.db
#   select batch_name, state, count(*) from ntuplizer_jobs group by state, batch_name;
#   select nretry, count(nretry) from ntuplizer_jobs group by nretry;
#   select log_file from ntuplizer_jobs WHERE state=="RUNNING" AND batch_name LIKE "TagAndProbeTrigger_220217_ParkingBPH4_Run2018D-05May2019promptD-v1_MINIAOD";
# To manage jobs you can use: -constraint 'JobBatchName == ""'
# condor_q -allusers -constraint 'regexp(".*B2DstMu.*",JobBatchName)'
# condor_hold -constraint 'regexp(".*SSDstMu_220301_ParkingBPH1.*",JobBatchName) && JobStatus == 2'
# Job status: 1=idle, 2=run, 5=hold
