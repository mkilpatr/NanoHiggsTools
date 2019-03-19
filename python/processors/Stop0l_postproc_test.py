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
from PhysicsTools.NanoSUSYTools.modules.TauMVAObjectsProducer import *
from PhysicsTools.NanoSUSYTools.modules.tauMVAProducer import *

mods = [
    #eleMiniCutID(),
    #Stop0lObjectsProducer("2017"),
    #DeepTopProducer("2017"),
    #Stop0lBaselineProducer("2017", isData=False, isFastSim=False),
    #UpdateGenWeight(args.crossSection, args.nEvents),
    #LLObjectsProducer(),
    #TauMVAObjectsProducer(),
    tauMVAProducer(),
]

#files=["/uscms/home/mkilpatr/nobackup/CMSSW_9_4_10/src/AnalysisMethods/macros/run/plots_19_01_30_smear/prod2017MC_NANO_Skim_original.root"]
#files=["/eos/uscms/store/user/lpcsusyhad/Stop_production/Summer16_80X_v2_NanAOD_MC/PostProcess_v1/TTbar_HT-600to800/TTbar_HT-600to800_0.root"]
#files=["/uscms_data/d3/lpcsusyhad/benwu/Moriond2019/TestNanoAOD/CMSSW_10_4_X_2018-12-11-2300/src/prod2017MC_NANO.root"]
#files=["root://cmseos.fnal.gov//store/user/benwu/Stop18/NtupleSyncMiniAOD/NanoSUSY/2018Xmas/prod2017MC_NANO.root"]
#files=["/eos/uscms/store/user/mkilpatr/13TeV/tauMVA/prod2017MC_NANO_Skim.root"]
files=["root://cmseos.fnal.gov//store/user/mkilpatr/13TeV/tauMVA/TTbarSingleLepT_2017_0.root"]
#files = []
#lines = open("samples_single.txt").readlines()
#for line in lines:
#    files.append(line.strip())
p=PostProcessor(".",files,cut=None, branchsel=None, outputbranchsel="keep_and_drop.txt", modules=mods,provenance=False)
#p=PostProcessor(".",files,cut=None, branchsel=None, outputbranchsel="keep_and_drop_QCD.txt", outputbranchselsmear="keep_and_drop_tauMVA.txt",modules=mods,provenance=False)
#p=PostProcessor(".",files,cut=None, branchsel=None, outputbranchsel="keep_and_drop_tauMVA.txt", typeofprocess="tau",modules=mods,provenance=False)
p.run()
