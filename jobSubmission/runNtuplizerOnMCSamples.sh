# Usage: jobSubmission/runNtuplizerOnMCSamples.sh

declare -a arr=(
"CP_General_MuEnriched_HardQCDall_ptHat15to20_TuneCP5_13TeV-pythia8"
"CP_General_MuEnriched_HardQCDall_ptHat20to30_TuneCP5_13TeV-pythia8"
"CP_General_MuEnriched_HardQCDall_ptHat30to50_TuneCP5_13TeV-pythia8"
)

for i in "${arr[@]}"
do
  process="$i"
  echo $process
  python jobSubmission/submitCMSSWCondorJobs.py -i production/inputFiles_$process.txt -o /storage/user/ocerri/BPhysics/data/cmsMC_private/$process/ntuples_B2DstMu/out_CAND.root -c config/cmssw_centralMC_Tag_B_MuDst-PiPiK_106x.py --maxtime 30m -N 20 -f

  echo "-----------------------"
  echo " "
  sleep 1
done
