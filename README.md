# NanoHiggs-tools
Postprocessing script for Stop 0L analysis

See the [readme](python/processors/Condor/README.md) in "NanoHiggsTools/python/processors/Condor" for specific instructions for condor submission.

### Set up CMSSW

```tcsh
source /cvmfs/cms.cern.ch/cmsset_default.csh
setenv SCRAM_ARCH slc6_amd64_gcc700
cmsrel CMSSW_10_2_9
cd CMSSW_10_2_9/src/
cmsenv
```

### Set up NanoHiggsTools framework
```tcsh
cd $CMSSW_BASE/src
cmsenv
git clone https://github.com/cms-tau-pog/TauIDSFs TauPOG/TauIDSFs
git clone -b Stop0l git@github.com:susy2015/nanoAOD-tools.git PhysicsTools/NanoAODTools
# For condor submission check the specific tag checkout instructions in [readme](python/processors/Condor/README.md)
git clone git@github.com:mkilpatr/NanoHiggsTools.git PhysicsTools/NanoHiggsTools
scram b
```

### Run local MC test
```tcsh
python Stop0l_postproc.py -i file:[input file] -s [MC sample name (from sampleSet cfg file)] -e [year]
```

### Process the MC Tuples on Condor to create higher level variables
```
python SubmitLPC.py -f ../Stop0l_postproc_diHiggs.py -c ../sampleSets_PreProcessed_2018.cfg -o /store/user/USER/outputdir/ -e 2018
```

### Run MC post-processing to use different selection methods
```
python SubmitLPC.py -f ../Stop0l_postproc_diHiggs.py -c ../sampleSets_PostProcessed_2018.cfg -o /store/user/mUSER/outputdir/ -e 2018 -p [dihiggs, cut, json, lund] -r [emu, muhad, ehad, hadhad, ""]
```

### Submit all regions for JSON generation and ParticleNet/LundNet training
```
cd Condor/
. ./submitRegion.sh <output director on EOS>
```
