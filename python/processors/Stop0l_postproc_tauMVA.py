#!/usr/bin/env python
import os, sys
import argparse
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
from importlib import import_module
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoSUSYTools.modules.eleMiniCutIDProducer import *
from PhysicsTools.NanoSUSYTools.modules.Stop0lObjectsProducer import *
from PhysicsTools.NanoSUSYTools.modules.Stop0lBaselineProducer import *
from PhysicsTools.NanoSUSYTools.modules.DeepTopProducer import *
from PhysicsTools.NanoSUSYTools.modules.updateEvtWeight import *
from PhysicsTools.NanoSUSYTools.modules.lepSFProducer import *
from PhysicsTools.NanoSUSYTools.modules.updateJetIDProducer import *
from PhysicsTools.NanoAODTools.postprocessing.modules.common.puWeightProducer import *
from PhysicsTools.NanoSUSYTools.modules.tauMVAProducer import *
from PhysicsTools.NanoSUSYTools.modules.PDFUncertaintyProducer import *
from PhysicsTools.NanoSUSYTools.modules.GenPartFilter import GenPartFilter
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jecUncertainties import jecUncertProducer

DataDepInputs = {
    "2016" : { "pileup": "Cert271036_284044_23Sep2016ReReco_Collisions16.root"
   },
    "2017" : { "pileup": "Cert294927_306462_EOY2017ReReco_Collisions17.root"
   },
    "2018" : { "pileup": "Cert314472_325175_PromptReco_Collisions18.root"
   }
}

def main(args):
    isdata = False
    isfastsim = False
    if "False" in args.isData:
        isdata = False
    else:
        isdata = True
    if "False" in args.isFastSim:
        isfastsim = False
    else:
        isfastsim = True

    mods = [
        eleMiniCutID(),
        Stop0lObjectsProducer(args.era),
        DeepTopProducer(args.era),
        Stop0lBaselineProducer(args.era, isData=isdata, isFastSim=isfastsim),
        UpdateEvtWeight(isdata, args.crossSection, args.nEvents),
	tauMVAProducer(),
    ]
    if args.era == "2018":
        mods.append(UpdateJetID(args.era))

    #~~~~~ For MC ~~~~~
    if not isdata:
        pufile = "%s/src/PhysicsTools/NanoSUSYTools/data/pileup/%s" % (os.environ['CMSSW_BASE'], DataDepInputs[args.era]["pileup"])
        mods += [
            # jecUncertProducer(DataDepInputs[args.era]["JECU"]),
            #PDFUncertiantyProducer(isdata),
            # lepSFProducer(args.era),
            puWeightProducer("auto", pufile, "pu_mc","pileup", verbose=False),
            # statusFlag 0x2100 corresponds to "isLastCopy and fromHardProcess"
            # statusFlag 0x2080 corresponds to "IsLastCopy and isHardProcess"
            GenPartFilter(statusFlags = [0x2100, 0x2080]),
        ]

    #files = ["/eos/uscms/store/user/lpcsusyhad/Stop_production/Summer16_94X_v3/PreProcessed_22Feb2019/TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/2016_MINIAODv3_RunIISummer16MiniAODv3-PUMoriond17_94X_v3-v2-ext1/190225_171125/0000/prod2016MC_NANO_1-1.root"]
    #files = ["/eos/uscms/store/group/lpcsusyhad/Stop_production/Autumn18_102X_v1/PreProcessed_22March2019//ZJetsToNuNu_HT-100To200_13TeV-madgraph/2018_Data_RunIIAutumn18MiniAOD-102X_v15-v1/190325_045752/0000/prod2018MC_NANO_1-1.root"]
    files = []
    lines = open(args.inputfile).readlines()
    for line in lines:
        files.append(line.strip())


    p=PostProcessor(args.outputfile,files,cut="MET_pt > 100 & nJet >= 2", branchsel=None, outputbranchsel="keep_and_drop_LL.txt", modules=mods,provenance=False)
    #p=PostProcessor(args.outputfile,files,cut="MET_pt > 100 & nJet >= 2", branchsel=None, outputbranchsel="keep_and_drop.txt", modules=mods,provenance=False)
    #p=PostProcessor(args.outputfile,files,cut=None, branchsel=None, outputbranchsel="keep_and_drop.txt", modules=mods,provenance=False)
    p.run()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='NanoAOD postprocessing.')
    parser.add_argument('-i', '--inputfile',
        default = "testing.txt",
        help = 'Path to the input filelist.')
    parser.add_argument('-o', '--outputfile',
                        default="./",
                        help = 'Path to the output file location.')
    parser.add_argument('-e', '--era',
        default = "2016", help = 'Year of production')
    parser.add_argument('-f', '--isFastSim', default = False)
    parser.add_argument('-d', '--isData', default = False)
    parser.add_argument('-c', '--crossSection',
                        type=float,
                        default = 1,
                        help = 'Cross Section of MC')
    parser.add_argument('-n', '--nEvents',
                        type=float,
                        default = 1,
                        help = 'Number of Events')
    args = parser.parse_args()
    main(args)
