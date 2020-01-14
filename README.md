# NanoSUSY-tools
Postprocessing script for Stop 0L analysis

See the [readme](python/processors/Condor/README.md) in "NanoSUSYTools/python/processors/Condor" for specific instructions for condor submission.

### Set up CMSSW

```tcsh
source /cvmfs/cms.cern.ch/cmsset_default.csh
setenv SCRAM_ARCH slc6_amd64_gcc700
cmsrel CMSSW_10_2_9
cd CMSSW_10_2_9/src/
cmsenv
```

### Set up NanoSusyTools framework
```tcsh
cd $CMSSW_BASE/src
cmsenv
git clone https://github.com/cms-tau-pog/TauIDSFs TauPOG/TauIDSFs
git clone -b Stop0l git@github.com:susy2015/nanoAOD-tools.git PhysicsTools/NanoAODTools
# For condor submission check the specific tag checkout instructions in [readme](python/processors/Condor/README.md)
git clone -b dev_v5 git@github.com:susy2015/NanoSUSY-tools.git PhysicsTools/NanoSUSYTools
git clone -b Stop0l_NanoAOD_production_V3.1 git@github.com:susy2015/TopTagger.git
scram b
cd $CMSSW_BASE/src/TopTagger/TopTagger/test
./configure
make
cmsenv
cd $CMSSW_BASE/src/PhysicsTools/NanoSUSYTools/python/processors
getTaggerCfg.sh -n -t DeepResolved_DeepCSV_GR_nanoAOD_2016_v1.0.3
getTaggerCfg.sh -n -t DeepResolved_DeepCSV_GR_nanoAOD_2017_v1.0.3
getTaggerCfg.sh -n -t DeepResolved_DeepCSV_GR_nanoAOD_2018_v1.0.3
```

### Run local MC test
```tcsh
python Stop0l_postproc.py -i file:[input file] -s [MC sample name (from sampleSet cfg file)] -e [year]
```

### Run local Data test
```tcsh
python Stop0l_postproc.py -i file:[input file] -d [data period] -e [year]
```


## To Do:
* PDF uncertainty module (Done)
    * weights stored in NanoAOD accordingly already, need code to extract the envelope
* lepton SF module follow [SUSY](https://twiki.cern.ch/twiki/bin/viewauth/CMS/SUSLeptonSF#Scale_Factors_for_SUSY_IDs) (Done with latest SF)
    * Need to update based upon the existing example code from NanoAOD-tools
* PU reweighting  (Done, to be tested)
    * Example code existed, but need recompute the pileup distribution from data
* btag SF update (instead of the btag weight stored during production) (done)
    * Need follow up with which method to be apply (iterative fit for DeepCSV?)
* JEC uncertainty update (done)
    * Need to understand the existing tool in NanoAOD-Tools
* DeepAK8/DeepResolved SF
    * Are they available yet?
* update Jet ID : [JetID Twiki](https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID13TeVRun2018) (Done with temp fix )
    * Do we need it in post-processing, or wait for next production?
* L1EcalPrefiring [twiki]* (https://twiki.cern.ch/twiki/bin/viewauth/CMS/L1ECALPrefiringWeightRecipe#Call_the_producer_in_your_config)
    * Rumor is it will be included in next NanoAOD
    * If not, we will make a module for it
* Various systematics 
* Trigger path and efficiency:
    * Once Hui's study is finalized, we will store the bit and efficiency + systematic

Smearing QCD Notes NANOAOD
In PhysicsTools/NanoSUSYTools/python/processors/
You need to create a QCD file with all of the qcd_orig files that you want to run over.
An example of all of these files is located in my area: /uscms_data/d3/{USER}/CMSSW_10_2_9/src/PhysicsTools/NanoSUSYTools/python/processors

> python Stop0l_postproc_QCD.py -p jetres

For Condor:
> python SubmitLPC_QCD.py -f ../Stop0l_postproc_QCD.py -c ../../../../../StopCfg/sampleSets_postProcess_2016_QCD.cfg -o /store/user/{USER}/13TeV/qcdsmearing_nanoaod -p jetres

hadd the output files together with the name "jetResSkim_combined_filtered_CHEF_NANO.root"

> mv jetResSkim_combined_filtered_CHEF_NANO.root /eos/uscms/store/user/{USER}/13TeV/qcdsmearing_nanoaod/.

You need to get the jet Res for the Tail Smear
> root -l -b -q ../rootlogon.C GetJetResForTailSmear.C+

> mkdir skims

> cd skims

> cp ../JetResDiagnostic.C .

Use JetResDiagnostic.C to create the pngs from the previous command ^
> root -l -b -q ../../rootlogon.C JetResDiagnostic.C+

> python Stop0l_postproc_QCD.py -p smear

For Condor:
> python SubmitLPC_QCD.py -f ../Stop0l_postproc_QCD.py -c ../../../../../StopCfg/sampleSets_PreProcessed_2016_QCD.cfg -o /store/user/{USER}/13TeV/qcdsmearing_nanoaod/ -i /store/user/{USER}/13TeV/qcdsmearing_nanoaod/resTailOut_combined_filtered_CHEF_puWeight_weight_WoH_NORMALIZED_NANO.root -p smear -m 4

After QCD smearing you need to run the add weight part of the NTuples.
You need to now create trees to run over and create the SF for the files. 
Once you create the trees from this file you need to have the met_tree.root, ttbarplusw_tree.root, and qcd_tree.root in the same directory.
You need to link the directory to the input files in the MakeQCDRespTailSF.C in the AnalysisMethods/macros/JetMETStudies/ directory.

> python Stop0l_postproc_QCD.py -p qcdsf

For Condor:
> python SubmitLPC_QCD.py -f ../Stop0l_postproc_QCD.py -c ../../../../../StopCfg/sampleSets_PostProcessed_2016_QCD_SF.cfg -o /store/user/mkilpatr/13TeV/nanoaod_QCDSF/ -e 2016 -p qcdsf

> cd AnalysisMethods/macros/JetMETStudies/

> root -l -b -q ../rootlogon.C MakeQCDRespTailSF_NANO.C+

This will run over the input files and calculate the correct SF for the QCD. Once it finishes there will be an output root file. You need to put this file in the data directory under corrections/ and in the right year. 

After getting the SF the QCD needs to be postprocessed one last time.
