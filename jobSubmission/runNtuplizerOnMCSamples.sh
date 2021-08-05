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

declare -a processes=(
    "CP_General_BdToJpsiKstar_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen"
    "BParking_Bd_JpsiKst_SoftQCDnonD_scale5_TuneCP5_HELAMP_PUc2_10-2-3"
    "CP_General_BuToJpsiK_BMuonFilter_TuneCP5_13TeV-pythia8-evtgen"
    "CP_BdToDstarMuNu_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen"
    "CP_BdToDstarTauNu_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen"
    "CP_BuToMuNuDstPi_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen"
    "CP_BdToMuNuDstPi_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen"
    "CP_BdToMuNuDstPiPi_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen"
    "CP_BuToMuNuDstPiPi_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen"
    "CP_BuToTauNuDstPi_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen"
    "CP_BdToTauNuDstPi_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen"
    "CP_BdToTauNuDstPiPi_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen"
    "CP_BuToTauNuDstPiPi_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen"
    "CP_BsToMuNuDstK_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen"
    "CP_BsToTauNuDstK_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen"
    "CP_BdToDstDu_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen"
    "CP_BdToDstDd_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen"
    "CP_BdToDstDs_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen"
    "CP_BuToDstDu_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen"
    "CP_BuToDstDd_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen"
    "CP_BsToDstDs_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen"
    "BP_Tag-Probe_B0_JpsiKst_Hardbbbar_evtgen_HELAMP_PUc0_10-2-3"
    "BP_Tag_B0_MuNuDmst_Hardbbbar_evtgen_ISGW2_PUc0_10-2-3"
)

output=$HOME/BPhysics/data/cmsMC

for process in "${processes[@]}"; do
    echo $process
    output_dir=$output/$process/ntuples/B2DstMu
    mkdir -p $output_dir
    python jobSubmission/create-condor-jobs -i production/inputFiles_$process.txt -o $output_dir/out_CAND.root -c config/cmssw_centralMC_Tag_B_MuDst-PiPiK.py -N 5
    sleep 1
done
