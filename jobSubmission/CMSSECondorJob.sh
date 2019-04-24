#!/bin/bash

ntuplizer_loc=$1
config=$2
inputFiles=$3
outputFile=$4


cd $ntuplizer_loc
eval `scramv1 runtime -sh`
cmsRun $config inputFiles=$inputFiles outputFile=$outputFile
