files=''
outfile = ''

#load specific moduels
from FWLite_plugins.BPHTriggerPath import BPHTriggerPath

exe_seq = [
           BPHTriggerPath(mu_charge=1, filter=True, produce_output=False)
          ]
