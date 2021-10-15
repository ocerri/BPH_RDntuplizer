# files=['/afs/cern.ch/user/o/ocerri/cernbox/BPhysics/data/cmsMC_private/BPH_Tag-Bm_D0kpmunu_Probe-Bp_D0kpmunu_NoPU_10-2-3_v0/BPH_Tag-Bm_D0kpmunu_Probe-Bp_D0kpmunu_MINIAODSIM.root']
#
# outfile=files[0].replace('MINIAODSIM.root', 'BPHMC.root')

from glob import glob
files=glob('/afs/cern.ch/user/o/ocerri/cernbox/BPhysics/data/cmsMC_private/BPH_Tag-Bm_D0kpmunu_Probe-Bp_D0kptaunu_tau2mununu_NoPU_10-2-3_v0/jobs_out/*_MINIAODSIM_*.root')

outfile=files[0].replace('jobs_out/', '')
outfile=outfile.replace('MINIAODSIM_1.root', 'BPHMC_merged.root')


#load specific moduels
from FWLite_plugins.MC_Process import MC_Process
from FWLite_plugins.BPHTriggerPath import BPHTriggerPath
from FWLite_plugins.RECO2MC_matching import RECO2MC_matching
from FWLite_plugins.EasyRecoQuantities import EasyRecoQuantities

exe_seq = [
           BPHTriggerPath(mu_charge=-1, filter=True),
           MC_Process(probe_pdgId=521, tag_pdgId=-521),
           RECO2MC_matching(verbose=False),
           EasyRecoQuantities(meas='RD', verbose=False)
          ]
