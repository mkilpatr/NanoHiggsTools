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

DataDepInputs = {
    "2016" : { "pileup": "Cert271036_284044_23Sep2016ReReco_Collisions16.root"
   },
    "2017" : { "pileup": "Cert294927_306462_EOY2017ReReco_Collisions17.root"
   },
    "2018" : { "pileup": "Cert314472_325175_PromptReco_Collisions18.root"
   }
}

def main(args):
    # isdata = False
    # isfastsim = False
    if "False" in args.isData:
        isdata = False
    else:
        isdata = True
    if "False" in args.isFastSim:
        isfastsim = False
    else:
        isfastsim = True

    mods = [
	Stop0lObjectsProducer(args.era),
	Stop0lBaselineProducer(args.era, isData=isdata, isFastSim=isfastsim),
	LLObjectsProducer(args.era),
    ]

    files = []
    if len(args.inputfile) > 5 and args.inputfile[0:5] == "file:":
        #This is just a single test input file
        files.append(args.inputfile[5:])
    else:
        #this is a file list
        with open(args.inputfile) as f:
            files = [line.strip() for line in f]

    #p=PostProcessor(args.outputfile,files,cut=None, branchsel=None, outputbranchsel="keep_and_drop_tauMVA.txt", typeofprocess="tau", modules=mods,provenance=False)
    p=PostProcessor(args.outputfile,files,cut="MET_pt > 200 & nJet >= 2", branchsel=None, outputbranchsel="keep_and_drop_LL.txt", modules=mods,provenance=False)
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
        default = "2017", help = 'Year of production')
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
