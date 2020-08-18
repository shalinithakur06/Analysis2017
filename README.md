# Analysis
   
### Download and compile the package  ###  
* source /cvmfs/cms.cern.ch/cmsset_default.sh
* cmsrel CMSSW_10_4_0
* cd CMSSW_10_4_0/src
* cmsenv
* https://github.com/ravindkv/ExLep2016Analysis.git
* cd ExLep2016Analysis/src
* make clean 
* make
* cd .. 

### Compile and run, in one go ### 
* root -l 'runMe.C("Analyzer")'

### Submit condor batch jobs  ###

* cd condor
* ./RunCond_TIFR.sh ntupleT2Paths.txt
