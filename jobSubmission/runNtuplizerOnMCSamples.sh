# Usage: jobSubmission/runNtuplizerOnMCSamples.sh

# process=BPH_Tag-Probe_B0_JpsiKst-mumuKpi-kp_13TeV-pythia8_SoftQCD_PTFilter5_0p0-evtgen_SVV_PU0_10-2-3
# process=BPH_Tag-Probe_B0_JpsiKst-mumuKpi-kp_13TeV-pythia8_Hardbbbar_PTFilter5_0p0-evtgen_SVV_PU0_10-2-3
# process=BPH_Tag-Probe_B0_JpsiKst-mumuKpi-kp_13TeV-pythia8_Hardbbbar_PTFilter5_0p0-evtgen_SVV_PU20_10-2-3
# process=BPH_Tag-Probe_B0_JpsiKst-mumuKpi-kp_13TeV-pythia8_Hardbbbar_PTFilter5_0p0-evtgenFSR_SVV_PU20_10-2-3
# process=BP_Tag-Probe_B0_JpsiKst_Hardbbbar_evtgen_HELAMP_PUc0_10-2-3
#
# python jobSubmission/submitCMSSWCondorJobs.py -i production/inputFiles_$process.txt -o /storage/user/ocerri/BPhysics/data/cmsMC_private/$process/ntuples_B2JpsiKst/out_CAND.root -c config/cmssw_privateMC_Tag_Mu-Probe-B0_JpsiKst-mumuKpi.py --maxtime 20m -N 50 -f
#
# python jobSubmission/submitCMSSWCondorJobs.py -i production/inputFiles_$process.txt -o /storage/user/ocerri/BPhysics/data/cmsMC_private/$process/ntuples_TagAndProbe/out_CAND.root -c config/cmssw_privateMC_TagAndProbe.py --maxtime 30m -N 60 -f

declare -a arr=(
"BP_Tag_B0_MuNuDmst_Hardbbbar_evtgen_ISGW2_PUc0_10-2-3"
"BP_Tag_B0_TauNuDmst_Hardbbbar_evtgen_ISGW2_PUc0_10-2-3"

##### "BP_Tag_B0_DmstHc_Hardbbbar_evtgen_ISGW2_PUc0_10-2-3"

"BP_Tag_B0_DstmDsp_Hardbbbar_evtgen_ISGW2_PUc0_10-2-3"
"BP_Tag_B0_DstmDp_Hardbbbar_evtgen_ISGW2_PUc0_10-2-3"
"BP_Tag_B0_DstmD0_Hardbbbar_evtgen_ISGW2_PUc0_10-2-3"

"BP_Tag_Bp_DstmHc_Hardbbbar_evtgen_ISGW2_PUc0_10-2-3"
"BP_Tag_Bm_DstmHc_Hardbbbar_evtgen_ISGW2_PUc0_10-2-3"
"BP_Tag_antiB0_DstmHc_Hardbbbar_evtgen_ISGW2_PUc0_10-2-3"
#
#
# "BP_Tag_Bp_MuNuDstst_Hardbbbar_evtgen_ISGW2_PUc0_10-2-3"
#
"BP_Tag_B0_MuNuDstst_Pi0_Hardbbbar_evtgen_ISGW2_PUc0_10-2-3"
"BP_Tag_B0_DmstPi0MuNu_Hardbbbar_evtgen_GR_PUc0_10-2-3"

"BP_Tag_Bp_MuNuDstst_PipPi0_Hardbbbar_evtgen_ISGW2_PUc0_10-2-3"
"BP_Tag_Bp_MuNuDstPipPi0_Hardbbbar_evtgen_PHSP_PUc0_10-2-3"

"BP_Tag_B0_MuNuDstst_PipPim_Hardbbbar_evtgen_ISGW2_PUc0_10-2-3"
"BP_Tag_B0_MuNuDstPipPim_Hardbbbar_evtgen_PHSP_PUc0_10-2-3"

"BP_Tag_B0_MuNuDstst_Pi0Pi0_Hardbbbar_evtgen_ISGW2_PUc0_10-2-3"
#
# "BP_Tag_B0_MuNuDstPiPiPi_Hardbbbar_evtgen_PHSP_PUc0_10-2-3"

# "BP_Tag_Bp_MuNuDstPiPiPi_Hardbbbar_evtgen_PHSP_PUc0_10-2-3"
#
#
# "BP_Tag_Bp_TauNuDstst_Pip_Hardbbbar_evtgen_ISGW2_PUc0_10-2-3"
# "BP_Tag_B0_TauNuDstst_Pi0_Hardbbbar_evtgen_ISGW2_PUc0_10-2-3"
)

for i in "${arr[@]}"
do
  process="$i"
  echo $process
  python jobSubmission/submitCMSSWCondorJobs.py -i production/inputFiles_$process.txt -o /storage/user/ocerri/BPhysics/data/cmsMC_private/$process/ntuples_B2DstMu/out_CAND.root -c config/cmssw_privateMC_Tag_B0_MuDmst-pD0bar-kp.py --maxtime 30m -N 100 -f

##### python jobSubmission/submitCMSSWCondorJobs.py -i /storage/user/ocerri/BPhysics/data/cmsMC_private/$process/jobs_out/out_MINIAODSIM_*.root -o /storage/user/ocerri/BPhysics/data/cmsMC_private/$process/jobs_B2DstMu/out_CAND.root -c config/cmssw_privateMC_Tag_B0_MuDmst-pD0bar-kp.py --maxtime 20m -N 20 -f ####

# python jobSubmission/submitCMSSWCondorJobs.py -i production/inputFiles_$process.txt -o /storage/user/ocerri/BPhysics/data/cmsMC_private/$process/ntuples_TagAndProbe_Bp_MuNuDstst/out_CAND.root -c config/cmssw_privateMC_TagAndProbe_Bp_MuNuDstst.py --maxtime 30m -N 100 -f
  echo "-----------------------"
  echo " "
  sleep 1
done
