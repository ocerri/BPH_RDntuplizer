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
    "CP_General_BdToJpsiKstar_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen"
    "CP_General_BuToJpsiK_BMuonFilter_TuneCP5_13TeV-pythia8-evtgen"
)

output=$HOME/BPhysics/data/cmsMC

for process in "${processes[@]}"; do
    echo $process
    output_dir=$output/$process/ntuples/B2DstMu
    if [ -d "$output_dir" ]; then
	echo "skipping $process because $output_dir exists";
    else
	mkdir -p $output_dir
	python jobSubmission/create-condor-jobs -i production/inputFiles_$process.txt -o $output_dir/out_CAND.root -c config/cmssw_centralMC_Tag_B_MuDst-PiPiK.py -N 5
	sleep 1
    fi
done
