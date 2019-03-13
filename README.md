# NanoSUSY-tools
Postprocessing script for Stop 0L analysis

### Set up CMSSW

```tcsh
source /cvmfs/cms.cern.ch/cmsset_default.csh
setenv SCRAM_ARCH slc6_amd64_gcc700
cmsrel CMSSW_10_2_6
cd CMSSW_10_2_6/src/
cmsenv
```

### Set up CMSSW
```tcsh
cd $CMSSW_BASE/src
cmsenv
git clone -b Stop0l git@github.com:susy2015/nanoAOD-tools.git PhysicsTools/NanoAODTools
git clone git@github.com:susy2015/NanoSUSY-tools.git PhysicsTools/NanoSUSYTools
scram b
```


## To Do:
* PDF uncertainty module 
    * weights stored in NanoAOD accordingly already, need code to extract the envelope
* lepton SF module follow [SUSY](https://twiki.cern.ch/twiki/bin/viewauth/CMS/SUSLeptonSF#Scale_Factors_for_SUSY_IDs)
    * Need to update based upon the existing example code from NanoAOD-tools
* PU reweighting 
    * Example code existed, but need recompute the pileup distribution from data
* btag SF update (instead of the btag weight stored during production)
    * Need follow up with which method to be apply (iterative fit for DeepCSV?)
* JEC uncertainty update
    * Need to understand the existing tool in NanoAOD-Tools
* DeepAK8/DeepResolved SF
    * Are they available yet?
* update Jet ID : [JetID Twiki](https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID13TeVRun2018)
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
An example of all of these files is located in my area: /uscms_data/d3/mkilpatr/CMSSW_10_2_9/src/PhysicsTools/NanoSUSYTools/python/processors

> python Stop0l_postproc_res.py

For Condor:
> python SubmitLPC.py -f ../Stop0l_postproc_res.py -c ../../../../../StopCfg/sampleSets_postProcess_2016_QCD.cfg -o /eos/uscms/store/user/{USER}/13TeV/qcdsmearing_nanoaod

hadd the output files together with the name "jetResSkim_combined_filtered_CHEF_NANO.root"

> mv jetResSkim_combined_filtered_CHEF_NANO.root /eos/uscms/store/user/{USER}/13TeV/qcdsmearing_nanoaod/.

You need to get the jet Res for the Tail Smear
> root -l -b -q ../rootlogon.C GetJetResForTailSmear.C+

> mkdir skims
> cd skims
> cp ../JetResDiagnostic.C .

Use JetResDiagnostic.C to create the pngs from the previous command ^
> root -l -b -q ../../rootlogon.C JetResDiagnostic.C+

> python Stop0l_postproc_QCD.py

For Condor:
Need 4 GB for the QCD smearing or else you run out of memory
> python SubmitLPC.py -f ../Stop0l_postproc_QCD.py -c ../../../../../StopCfg/sampleSets_preProcess_2016_QCD.cfg -o /eos/uscms/store/user/{USER}/13TeV/qcdsmearing_nanoaod/ -m 4

After QCD smearing you need to run the add weight part of the NTuples.
You need to now create trees to run over and create the SF for the files. 
Once you create the trees from this file you need to have the met_tree.root, ttbarplusw_tree.root, and qcd_tree.root in the same directory.
You need to link the directory to the input files in the MakeQCDRespTailSF.C in the AnalysisMethods/macros/JetMETStudies/ directory.

> cd AnalysisMethods/macros/JetMETStudies/
> root -l -b -q ../rootlogon.C MakeQCDRespTailSF.C+

This will run over the input files and calculate the correct SF for the QCD. Once it finishes there will be an output root file. You need to put this file in the data directory under corrections/ and in the right year. 

After getting the SF the QCD needs to be postprocessed one last time.
