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

ntuplesName=ntuples_B2DstMu_220209

config=config/cmssw_centralMC_Tag_Bd_MuDst-PiPiK.py

declare -a processes=(
    # "BP_Tag-Probe_B0_JpsiKst_Hardbbbar_evtgen_HELAMP_PUc0_10-2-3"
    # "BP_Tag_B0_MuNuDmst_Hardbbbar_evtgen_ISGW2_PUc0_10-2-3"
    #
    # "CP_General_BdToJpsiKstar_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen"
    # "CP_General_BuToJpsiK_BMuonFilter_TuneCP5_13TeV-pythia8-evtgen"
    # "BParking_Bd_JpsiKst_SoftQCDnonD_scale5_TuneCP5_HELAMP_PUc2_10-2-3"
    # "CP_General_MuEnriched_HardQCDall_TuneCP5_13TeV-pythia8"
    #
    # Central production --> Should be run N = 3
    "CP_BdToDstarMuNu_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen"
    "CP_BdToDstarTauNu_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen"
    "CP_BuToMuNuDstPi_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen"
    "CP_BdToMuNuDstPi_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen"
    "CP_BuToTauNuDstPi_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen"
    "CP_BdToTauNuDstPi_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen"
    "CP_BsToMuNuDstK_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen"
    "CP_BsToTauNuDstK_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen"
    "CP_BdToDstDu_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen"
    "CP_BdToDstDd_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen"
    "CP_BdToDstDs_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen"
    "CP_BuToDstDu_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen"
    "CP_BuToDstDd_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen"
    "CP_BsToDstDs_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen"
    #
    # Private production --> Should be run N = 100
    # "CP_BdToMuNuDstPiPi_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen_v3"
    # "CP_BuToMuNuDstPiPi_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen_v3"
    # "CP_BdToTauNuDstPiPi_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen"
    # "CP_BuToTauNuDstPiPi_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen"
)

for process in "${processes[@]}"; do
    echo $process
    output_dir=$outLoc/$process/$ntuplesName
    mkdir -p $output_dir
    python jobSubmission/create-condor-jobs -i production/inputFiles_$process.txt -o $output_dir/out_CAND.root -c $config -t $ntuplesName -N 3 --maxtime 120m
    sleep 1
done


# To actually submit the jobs run: ./jobSubmission/submit-condor-jobs --max-jobs 1000
# To check the status of the submission:
#   sqlite3 ~/state.db
#   sqlite3 -column -header ~/state.db
#   select batch_name, state, count(*) from ntuplizer_jobs group by state, batch_name;
#   select nretry, count(nretry) from ntuplizer_jobs group by nretry;
