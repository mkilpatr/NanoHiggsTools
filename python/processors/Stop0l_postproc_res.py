#!/usr/bin/env python
import os, sys
import ROOT
import argparse
ROOT.PyConfig.IgnoreCommandLineOptions = True
from importlib import import_module
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoSUSYTools.modules.JetResSkim import *

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
        JetResSkim(args.era),
    ]
    
    #files=["root://cmseos.fnal.gov//eos/uscms/store/user/lpcsusyhad/Stop_production/Summer16_80X_v2_NanAOD_MC//QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/2016_MC_RunIISummer16MiniAODv2-PUMoriond17_80X_v6-v1-ext1/190124_030402/0000/prod2016MC_NANO_1-33.root"]
    files = []
    lines = open(args.inputfile).readlines()
    for line in lines:
        files.append(line.strip())

    p=PostProcessor(args.outputfile,files,cut=None, branchsel=None, outputbranchsel="keep_and_drop_res.txt",typeofprocess="resp",modules=mods,provenance=False)
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
    parser.add_argument('-f', '--isFastSim', default = False,
                        help = "Input file is fastsim (Default: false)")
    parser.add_argument('-d', '--isData',    default = False,
                        help = "Input file is data (Default: false)")
    parser.add_argument('-s', '--sampleName', action="store",  type=str, default = "",
                        help = "Name of MC sample (from sampleSet file) (Default: )")
    parser.add_argument('-c', '--crossSection',
                        type=float,
                        default = 1,
                        help = 'Cross Section of MC to use for MC x-sec*lumi weight (Default: 1.0)')
    parser.add_argument('-n', '--nEvents',
                        type=float,
                        default = 1,
                        help = 'Number of events to use for MC x-sec*lumi weight (NOT the number of events to run over) (Default: 1.0)')
    parser.add_argument('-m', '--maxEvents',
                        type=int,
                        default = -1,
                        help = 'MAximum number of events to process (Default: all events)')
    args = parser.parse_args()
    main(args)
