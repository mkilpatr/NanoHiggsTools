#!/usr/bin/env python
import os, sys
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
from importlib import import_module
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoSUSYTools.modules.updateGenWeight import *
from PhysicsTools.NanoSUSYTools.modules.lepSFProducer import *
from PhysicsTools.NanoSUSYTools.modules.updateJetIDProducer import *
from PhysicsTools.NanoAODTools.postprocessing.modules.common.puWeightProducer import *
from PhysicsTools.NanoSUSYTools.modules.eleMiniCutIDProducer import *
from PhysicsTools.NanoSUSYTools.modules.Stop0lObjectsProducer import *
from PhysicsTools.NanoSUSYTools.modules.Stop0lBaselineProducer import *
from PhysicsTools.NanoSUSYTools.modules.DeepTopProducer import *
from PhysicsTools.NanoSUSYTools.modules.LLObjectsProducer import *

DataDepInputs = {
    "2016" : { "pileup": "Cert271036_284044_23Sep2016ReReco_Collisions16.root"
   },
    "2017" : { "pileup": "Cert294927_306462_EOY2017ReReco_Collisions17.root"
   },
    "2018" : { "pileup": "Cert314472_325175_PromptReco_Collisions18.root"
   }
}

isdata = False
isfastsim = False
#if "False" in args.isData:
#    isdata = False
#else:
#    isdata = True
#if "False" in args.isFastSim:
#    isfastsim = False
#else:
#    isfastsim = True

mods = [
    eleMiniCutID(),
    Stop0lObjectsProducer("2017"),
    DeepTopProducer("2017"),
    Stop0lBaselineProducer("2017", isData=isdata, isFastSim=isfastsim),
    #UpdateGenWeight(isdata, args.crossSection, args.nEvents)
    LLObjectsProducer(),
]
if "2017" == "2018":
    mods.append(UpdateJetID("2017"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ For MC ~~~~~
if not isdata:
    pufile = "%s/src/PhysicsTools/NanoSUSYTools/data/pileup/%s" % (os.environ['CMSSW_BASE'], DataDepInputs["2017"]["pileup"])
    mods += [
        lepSFProducer("2017"),
        puWeightProducer("auto", pufile, "pu_mc","pileup", verbose=False)
    ]


#files=["/uscms/home/mkilpatr/nobackup/CMSSW_9_4_10/src/AnalysisMethods/macros/run/plots_19_01_30_smear/prod2017MC_NANO_Skim_original.root"]
#files=["/uscms_data/d3/lpcsusyhad/benwu/Moriond2019/TestNanoAOD/CMSSW_10_4_X_2018-12-11-2300/src/prod2017MC_NANO.root"]
#files=["root://cmseos.fnal.gov//store/user/benwu/Stop18/NtupleSyncMiniAOD/NanoSUSY/2018Xmas/prod2017MC_NANO.root"]
#files=["/eos/uscms/store/user/mkilpatr/13TeV/tauMVA/prod2017MC_NANO_Skim.root"]
files=["/eos/uscms/store/user/lpcsusyhad/Stop_production/Fall17_94X_v2_NanAOD_MC/PreProcessed_15Jan2019/TTJets_SingleLeptFromT_TuneCP5_13TeV-madgraphMLM-pythia8/2017_MC_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_v14-v1/190111_191501/0000/prod2017MC_NANO_1.root"]
#p=PostProcessor(".",files,cut=None, branchsel=None, outputbranchsel="keep_and_drop_QCD.txt", outputbranchselsmear="keep_and_drop_tauMVA.txt",modules=mods,provenance=False)
#p=PostProcessor(".",files,cut=None, branchsel=None, outputbranchsel="keep_and_drop_tauMVA.txt", typeofprocess="tau",modules=mods,provenance=False)
p=PostProcessor(".",files,cut=None, branchsel=None, outputbranchsel="keep_and_drop.txt", modules=mods,provenance=False)
p.run()
