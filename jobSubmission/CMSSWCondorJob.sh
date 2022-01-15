#!/bin/bash

ntuplizer_loc=$1
config=$2
inputFiles=$3
outputFile=$4
fjrFile=$5


cd $ntuplizer_loc
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scramv1 runtime -sh`
cmsRun --jobreport=$fjrFile $config inputFiles=$inputFiles outputFile=$outputFile
exitcode=$?
echo "CMSSW exit code: $exitcode"
echo "========= DONE ========="
exit $exitcode
