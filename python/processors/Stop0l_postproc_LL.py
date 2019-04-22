#!/usr/bin/env python
import os, sys
import argparse
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
from importlib import import_module
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoSUSYTools.modules.Stop0lObjectsProducer import *
from PhysicsTools.NanoSUSYTools.modules.Stop0lBaselineProducer import *
from PhysicsTools.NanoSUSYTools.modules.LLObjectsProducer import *
from PhysicsTools.NanoSUSYTools.modules.tauMVA import *

DataDepInputs = {
    "2016" : { "pileup": "Cert271036_284044_23Sep2016ReReco_Collisions16.root",
               "JECU": "Summer16_07Aug2017_V11_MC"
               },
    "2017" : { "pileup": "Cert294927_306462_EOY2017ReReco_Collisions17.root",
               "JECU": "Fall17_17Nov2017_V32_MC"
               },
    "2018" : { "pileup": "Cert314472_325175_PromptReco_Collisions18.root",
                #The 2018 files is actually a softlink to this file
               "JECU": "Fall17_17Nov2017_V32_MC"
               }
}

def main(args):
    isdata = args.isData
    isfastsim = args.isFastSim
    print(isdata, isfastsim)

    if isdata and isfastsim:
        print "ERROR: It is impossible to have a dataset that is both data and fastsim"
        exit(0)

    if not args.era in DataDepInputs.keys():
        print "ERROR: Era \"" + args.era + "\" not recognized"
        exit(0)

    mods = [
	Stop0lObjectsProducer(args.era),
	#tauMVA(),
	Stop0lBaselineProducer(args.era, isData=isdata, isFastSim=isfastsim),
	LLObjectsProducer(args.era),
    ]

    files=["/eos/uscms/store/user/lpcsusyhad/Stop_production/Fall17_94X_v2_NanAOD_MC/PostProcessed_15Jan2019_v1/TTbar_HT-600to800_2017/TTbar_HT-600to800_2017_0.root"]
    #files=["root://cmseos.fnal.gov//eos/uscms/store/user/lpcsusyhad/Stop_production/Summer16_94X_v3/PostProcessed_22Feb2019_v2p2/TTbarDiLep_2016/TTbarDiLep_2016_10.root"]
    #files = []
    #if len(args.inputfile) > 5 and args.inputfile[0:5] == "file:":
    #    #This is just a single test input file
    #    files.append(args.inputfile[5:])
    #else:
    #    #this is a file list
    #    with open(args.inputfile) as f:
    #        files = [line.strip() for line in f]

    #p=PostProcessor(args.outputfile,files,cut=None, branchsel=None, outputbranchsel="keep_and_drop_tauMVA.txt", typeofprocess="tau", modules=mods,provenance=False)
    p=PostProcessor(args.outputfile,files,cut="MET_pt > 200 & nJet >= 2", branchsel=None, outputbranchsel="keep_and_drop_LL.txt", modules=mods,provenance=False)
    #p=PostProcessor(args.outputfile,files,cut=None, branchsel=None, outputbranchsel="keep_and_drop.txt", modules=mods,provenance=False)
    p.run()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='NanoAOD postprocessing.')
    parser.add_argument('-i', '--inputfile',
        default = "testing.txt",
        help = 'Path to the input filelist. To run with a single file instead of a file list prepend the filepath with \"file:\" (Default: testing.txt)')
    parser.add_argument('-o', '--outputfile',
                        default="./",
                        help = 'Path to the output file location. (Default: .)')
    parser.add_argument('-e', '--era',
        default = "2017", help = 'Year of production')
    parser.add_argument('-f', '--isFastSim', action="store_true",  default = False,
                        help = "Input file is fastsim (Default: false)")
    parser.add_argument('-d', '--isData',  default = False,
                        help = "Input file is data (Default: false)")
    parser.add_argument('-c', '--crossSection',
                        type=float,
                        default = 1,
                        help = 'Cross Section of MC to use for MC x-sec*lumi weight (Default: 1.0)')
    parser.add_argument('-n', '--nEvents',
                        type=float,
                        default = 1,
                        help = 'Number of events to use for MC x-sec*lumi weight (NOT the number of events to run over) (Default: 1.0)')
    args = parser.parse_args()
    main(args)
