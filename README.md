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
wget https://hammer.physics.lbl.gov/Hammer-1.2.1-Source.tar.gz
tar -xzf Hammer-1.2.1-Source.tar.gz
mkdir Hammer-build
cd Hammer-build

export PATH=/cvmfs/sft.cern.ch/lcg/contrib/CMake/3.14.2/Linux-x86_64/bin:${PATH}
export BOOST_ROOT=/cvmfs/sft.cern.ch/lcg/releases/Boost/1.66.0-f50b5/x86_64-centos7-gcc7-opt/
export HEPMC_DIR=/cvmfs/sft.cern.ch/lcg/external/HepMC/2.06.08/x86_64-slc6-gcc48-opt


cmake -DCMAKE_INSTALL_PREFIX=../Hammer-install -DENABLE_TESTS=ON -DWITH_ROOT=OFF -DWITH_EXAMPLES=OFF -DINSTALL_EXTERNAL_DEPENDENCIES=ON -DWITH_PYTHON=OFF -DBUILD_SHARED_LIBS=ON -DFORCE_YAMLCPP_INSTALL=ON -DCMAKE_CXX_FLAGS="-pthread" -DBOOST_ROOT=/cvmfs/sft.cern.ch/lcg/releases/Boost/1.64.0-0809c/x86_64-centos7-gcc7-opt/ ../Hammer-1.2.1-Source
```
Change ln 79 of `Hammer-1.2.1-Source/src/Amplitudes/AmplBDstarDPiLepNu.cc` from `const FourMomentum& pPion = daughters[4].momentum();` to `const FourMomentum& pPion = pDstarmes - pDmes;`. Then continue to compilation.
```
make -j24; make install -j24

cd $CMSSW_BASE/lib/slc7_amd64_gcc700
cp ../../src/ntuplizer/Hammer-install/lib64/*.so.* ./
cp ../../src/ntuplizer/Hammer-install/lib64/Hammer/*.so.* ./
cd $CMSSW_BASE/src/ntuplizer/Hammer-build

ctest -V

cd ..
rm -rf Hammer-1.2.1-Source.tar.gz Hammer-build
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

### Submitting jobs to condor

To submit ntuplizer jobs to condor there is a two step process. First, we run
the script create-condor-jobs which creates the submission files and adds
entries to an sqlite database located at ~/state.db. For example, to run the
ntuplizer on all the BdToDstarMuNu Monte Carlo, you would run:

```console
$ process="CP_BdToDstarMuNu_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen"
$ python jobSubmission/create-condor-jobs -i production/inputFiles_$process.txt -o ~/BPhysics/$process/ntuples/out_CAND.root -c config/cmssw_centralMC_Tag_Bd_MuDst-PiPiK.py -N 5
```

You can then browse the submissions in the sqlite database by running:

```console
$ sqlite3 ~/state.db
sqlite> select timestamp, id, state from ntuplizer_jobs;
```

Note that by default if you run create-condor-jobs twice, it will only add new
submissions for files which haven't been processed yet. If you want to force it
to reprocess new files (for example with a different configuration file), you
need to pass the `--force` argument to the script.

The second step is to actually submit the jobs to condor. To submit the jobs you can run:

```console
$ ./jobSubmission/submit-condor-jobs --max-jobs 1000
```

This script is meant to be run from a cron job, so you can also set up a cron job like:

(Tony settings)
```
HOME=/storage/af/user/[username]
PATH=/usr/bin

0 * * * * source $HOME/.bashrc; ~/RDstAnalysis/CMSSW_10_2_3/src/ntuplizer/BPH_RDntuplizer/jobSubmission/submit-condor-jobs --max-jobs 100 --loglevel notice --logfile $HOME/submit.log
```

(Olmo settings)
```
HOME=/storage/af/user/ocerri
PATH=/usr/bin

0 * * * * source $HOME/.bash_profile; source /cvmfs/cms.cern.ch/cmsset_default.sh; cd $HOME/work/CMSSW_10_2_3/src; eval `scramv1 runtime -sh`; cd ntuplizer/BPH_RDntuplizer; ./jobSubmission/submit-condor-jobs --max-jobs 1000 --loglevel notice --logfile $HOME/submit.log 2>&1 >> $HOME/crontab.log
```

It will loop over the entries in the database, and for any that are in the
'NEW' state it will submit the job. For any jobs which complete, it will check
the log file to look to see that it successfully returned 0, and that the ROOT
file is not corrupted. If both of these checks pass it sets the state to
'SUCCESSFUL' and if not to 'FAILED'.

You can then check on any failed jobs by looking for them in the database:

```console
$ sqlite3 ~/state.db
sqlite> select timestamp, id, state from ntuplizer_jobs where state = 'FAILED';
```

If it's some random issue, you can tell submit-condor-jobs to retry it by doing:

```console
sqlite> update ntuplizer_jobs set state = 'RETRY' where id = [id];
```

Then, the next time you run submit-condor-jobs it will submit the job again up
to a maximum of a certain number of retries.

### Real data ntuplization

Real data ntuplization is done through CRAB. The submission script is [jobSubmission/submitCMSSWCrabJobs.py](jobSubmission/submitCMSSWCrabJobs.py).
While jobs are running the status can be checked with [jobSubmission/checkCrabStatus.py](jobSubmission/checkCrabStatus.py).
Finally the output can be collected using [jobSubmission/collectCrabJobs.py](jobSubmission/collectCrabJobs.py).
