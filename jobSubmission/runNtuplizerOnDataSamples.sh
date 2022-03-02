#!/usr/bin/env bash
# Script to submit jobs to condor for running the ntuplizer. To run it just
# run:
#
#     # sh ./jobSubmission/runNtuplizerOnDataSamples.sh
#
# Make sure you have updated the input files by running:
#
#     # cd production
#     # python createDatasetFileList.py
#
# first.

outLoc=/storage/af/group/rdst_analysis/BPhysics/data/cmsRD
maxRunTime=48h
nice=1
######## Trigger tag and probe ntuples #############
# ntuplesName=TagAndProbeTrigger_220217
# config=config/cmssw_cmsRD2018_TagAndProbeTrigger.py
# nFilesPerJob=40 # About 3h per job
# Check how much time it really took for these jobs!

######## Bd -> Jpsi K* ntuples #############
# ntuplesName=Bd2JpsiKst_220228
# config=config/cmssw_cmsRD2018_Bd_JpsiKst-mumuKpi.py
# nFilesPerJob=30
# Check how much time it really took for these jobs!

######## Main R(D*) analysis ntuples #############
# ntuplesName=B2DstMu_220220
# config=config/cmssw_cmsRD2018_Tag_Bd_MuDst-PiPiK.py
# nFilesPerJob=20 # About 4h per job
# Check how much time it really took for these jobs!

ntuplesName=SSDstMu_220301
config=config/cmssw_cmsRD2018_Tag_SS_MuDst-PiPiK.py
nFilesPerJob=30 # About 4h per job


for iPart in {1..5}; do
    echo Part $iPart
    output_dir=$outLoc/ParkingBPH${iPart}/Run2018D-05May2019promptD-v1_RDntuplizer_$ntuplesName
    echo $output_dir
    mkdir -p $output_dir
    python jobSubmission/create-condor-jobs -i production/inputFiles_ParkingBPH${iPart}_Run2018D-05May2019promptD-v1_MINIAOD.txt -o $output_dir/out_CAND.root -c $config -t $ntuplesName -N $nFilesPerJob --maxtime $maxRunTime --nice $nice
    sleep 1
done


# To actually submit the jobs run: ./jobSubmission/submit-condor-jobs --max-jobs 1000
# To check the status of the submission:
#   sqlite3 ~/state.db
#   sqlite3 -column -header ~/state.db
#   select batch_name, state, count(*) from ntuplizer_jobs group by state, batch_name;
#   select nretry, count(nretry) from ntuplizer_jobs group by nretry;
