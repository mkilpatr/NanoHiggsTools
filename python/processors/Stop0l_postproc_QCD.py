#!/usr/bin/env python
import os, sys
import ROOT
import argparse
ROOT.PyConfig.IgnoreCommandLineOptions = True
from importlib import import_module
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoSUSYTools.modules.qcdSmearProducer import *

#def main(args):
#    # isdata = False
#    # isfastsim = False
#    if "False" in args.isData:
#        isdata = False
#    else:
#        isdata = True
#    if "False" in args.isFastSim:
#        isfastsim = False
#    else:
#        isfastsim = True

mods = [
    qcdSmearProducer(),
]

files=["/eos/uscms/store/user/lpcsusyhad/Stop_production/Summer16_80X_v2_NanAOD_MC/PostProcess_v1/QCD_HT100to200/QCD_HT100to200_0.root"]
#files = []
#lines = open(args.inputfile).readlines()
#for line in lines:
#    files.append(line.strip())

#p=PostProcessor(args.outputfile,files,cut=None, branchsel=None, outputbranchsel="keep_and_drop_QCD.txt", outputbranchselsmear="keep_and_drop_smear.txt",typeofprocess="smear",modules=mods,provenance=False)
p=PostProcessor(".",files,cut=None, branchsel=None, outputbranchsel="keep_and_drop_QCD.txt", outputbranchselsmear="keep_and_drop_smear.txt",typeofprocess="smear",modules=mods,provenance=False)
p.run()

#if __name__ == "__main__":
#    parser = argparse.ArgumentParser(description='NanoAOD postprocessing.')
#    parser.add_argument('-i', '--inputfile',
#        default = "testing.txt",
#        help = 'Path to the input filelist.')
#    parser.add_argument('-o', '--outputfile',
#                        default="./",
#                        help = 'Path to the output file location.')
#    parser.add_argument('-e', '--era',
#        default = "2017", help = 'Year of production')
#    parser.add_argument('-f', '--isFastSim', default = False)
#    parser.add_argument('-d', '--isData', default = False)
#    parser.add_argument('-c', '--crossSection',
#                        type=float,
#                        default = 1,
#                        help = 'Cross Section of MC')
#    parser.add_argument('-n', '--nEvents',
#                        type=float,
#                        default = 1,
#                        help = 'Number of Events')
#    args = parser.parse_args()
#    main(args)
