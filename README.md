# BPH_RDntuplizer

Code to produce ntuples for the R(D*) analysis. They can be produced with FWLite or with the full CMSSW.

## Installing instructions
Has to run on lxplus.

```
cmsrel CMSSW_10_2_3
cd CMSSW_10_2_3/src
cmsenv
mkdir ntuplizer
cd ntuplizer
git clone https://github.com/ocerri/BPH_RDntuplizer.git
cd BPH_RDntuplizer
scram b -j12
```

## Running instructions

Running cmssw
```
cmsRun config/cmssw_myconfig.py
```

Running FWLite
```
python FWLite_master.py config/cfg_myconfig.py
```
