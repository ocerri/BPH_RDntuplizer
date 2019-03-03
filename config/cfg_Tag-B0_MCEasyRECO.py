from glob import glob
import re

files=''
outfile = ''

def makeOutname(l):
    if '/jobs_out/' in l:
        l = l.replace('/jobs_out/', '/')
    l = l.replace('MINIAODSIM', 'BPHMCEasyRECO')
    if not 'merge' in l:
        out = re.search('_[0-9]+\.root]', l)
        if hasattr(out, 'group'):
            l = l.replace(out.group(0), '.root')
    return l



#load specific moduels
from FWLite_plugins.MC_ProcessTag import MC_ProcessTag
from FWLite_plugins.BPHTriggerPath import BPHTriggerPath
from FWLite_plugins.RECO2MC_matching import RECO2MC_matching
from FWLite_plugins.EasyRecoQuantities import EasyRecoQuantities

exe_seq = [
           BPHTriggerPath(mu_charge=1, filter=True),
           MC_ProcessTag(tag_pdgId=511),
           RECO2MC_matching(verbose=False),
           EasyRecoQuantities(verbose=False)
          ]
