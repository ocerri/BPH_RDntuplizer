#!/bin/bash

ntuplizer_loc=$1
config=$2
inputFiles=$3
outputFile=$4


# cd $ntuplizer_loc
cd /storage/af/user/ocerri/work/CMSSW_10_6_25/src/ntuplizer/BPH_RDntuplizer # temporary fix
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scramv1 runtime -sh`
cmsRun $config inputFiles=$inputFiles outputFile=$outputFile
exitcode=$?
echo "CMSSW exit code: $exitcode"
echo "========= DONE ========="
exit $exitcode
