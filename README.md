# BPH_RDntuplizer

Code to produce ntuples for the R(D*) analysis. They can be produced with FWLite or with the full CMSSW.

## Installing instructions
Has to run on a CMS cluster (caltech T2, lxplus, etc.).

Create the CMSSW environment.
```
cmsrel CMSSW_10_2_3
cd CMSSW_10_2_3/src
cmsenv
mkdir ntuplizer
cd ntuplizer
```

Install [Hammer](https://gitlab.com/mpapucci/Hammer/-/tree/master)
```
export PATH=/cvmfs/sft.cern.ch/lcg/contrib/CMake/3.14.2/Linux-x86_64/bin:${PATH}
export BOOST_ROOT=/cvmfs/sft.cern.ch/lcg/releases/Boost/1.66.0-f50b5/x86_64-centos7-gcc7-opt/
export HEPMC_DIR=/cvmfs/sft.cern.ch/lcg/external/HepMC/2.06.08/x86_64-slc6-gcc48-opt

cp -r /storage/af/user/ocerri/public_html/Hammer-0.9.0-Source ./
mkdir Hammer-0.9.0-build
cd Hammer-0.9.0-build

cmake -DCMAKE_INSTALL_PREFIX=../Hammer-0.9.0 -DENABLE_TESTS=ON -DWITH_ROOT=OFF -DWITH_EXAMPLES=OFF -DINSTALL_EXTERNAL_DEPENDENCIES=ON -DWITH_PYTHON=OFF -DBUILD_SHARED_LIBS=ON ../Hammer-0.9.0-Source
make -j24
make install -j24

cd $CMSSW_BASE/lib/slc7_amd64_gcc700
cp ../../src/ntuplizer/Hammer-0.9.0/lib64/*.so.* ./
cp ../../src/ntuplizer/Hammer-0.9.0/lib64/Hammer/*.so.* ./
cd $CMSSW_BASE/src/ntuplizer/Hammer-0.9.0-build
ctest

cd ..
rm -rf Hammer-0.9.0-*
```


Fetch the ntuplizer and comple it all
```
git clone https://github.com/ocerri/BPH_RDntuplizer.git
cd BPH_RDntuplizer
scram b -j12
```

## Running instructions

Create list of files for each sample

```
cd production/
python createDatasetFileList.py
```


Running cmssw
```
cmsRun config/cmssw_myconfig.py
```

Running FWLite
```
python FWLite_master.py config/cfg_myconfig.py
```

To run the MC ntuplization jobs check [notes.txt](jobSubmission/notes.txt) file.

### Real data ntuplization

Real data ntuplization is done through CRAB. The submission script is [jobSubmission/submitCMSSWCrabJobs.py](jobSubmission/submitCMSSWCrabJobs.py).
While jobs are running the status can be checked with [jobSubmission/checkCrabStatus.py](jobSubmission/checkCrabStatus.py).
Finally the output can be collected using [jobSubmission/collectCrabJobs.py](jobSubmission/collectCrabJobs.py).
