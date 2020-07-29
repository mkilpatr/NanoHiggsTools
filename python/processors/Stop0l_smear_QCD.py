#!/usr/bin/env python
import os, sys
import ROOT
import argparse
ROOT.PyConfig.IgnoreCommandLineOptions = True
from importlib import import_module
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoSUSYTools.modules.qcdSmearProducer import *
from PhysicsTools.NanoSUSYTools.modules.JetResSkim import *
from PhysicsTools.NanoSUSYTools.modules.QCDObjectsProducer import *
from PhysicsTools.NanoSUSYTools.modules.qcdSFProducer import *
from PhysicsTools.NanoSUSYTools.modules.updateEvtWeightSmear import *
from PhysicsTools.NanoSUSYTools.modules.LLObjectsProducer import *

def main(args):
    isdata = len(args.dataEra) > 0
    isqcd = args.sampleName.startswith("QCD_")
    process = args.process

    mods = []
    if process == 'jetres':
	mods.append(JetResSkim(args.era))
    elif process == 'smear':
	mods.append(UpdateEvtWeightSmear(isdata, args.crossSection, args.nEvents, args.sampleName))
	mods.append(qcdSmearProducer())
    elif process == 'qcdsf':
	mods.append(QCDObjectsProducer(isQCD=isqcd, isData=isdata))
	mods.append(LLObjectsProducer(args.era, isData=isdata))
    elif process == 'sfcalc':
	mods.append(qcdSFProducer(args.era))
    
    files = []
    if len(args.inputfile) > 5 and args.inputfile[0:5] == "file:":
        #This is just a single test input file
        files.append(args.inputfile[5:])
    else:
        #this is a file list
        with open(args.inputfile) as f:
            files = [line.strip() for line in f]
    
    if process=='jetres':   p=PostProcessor(args.outputfile,files,cut=None, branchsel=None, outputbranchsel="keep_and_drop_res.txt",typeofprocess="resp",modules=mods,provenance=False,maxEvents=args.maxEvents)
    elif process=='smear':  p=PostProcessor(args.outputfile,files,cut=None, branchsel=None, outputbranchsel="keep_and_drop_QCD.txt", outputbranchselsmear="keep_and_drop_QCD.txt",typeofprocess="smear",modules=mods,provenance=False,maxEvents=args.maxEvents)
    #elif process=='qcdsf':  p=PostProcessor(args.outputfile,files,cut="Pass_MET & Pass_NJets20 & Pass_EventFilter & Pass_HT & Pass_JetID", branchsel=None, outputbranchsel="keep_and_drop.txt", modules=mods,provenance=False,maxEvents=args.maxEvents)
    elif process=='qcdsf':  p=PostProcessor(args.outputfile,files,cut="Pass_MET & nJet>=2", branchsel=None, outputbranchsel="keep_and_drop.txt", modules=mods,provenance=False,maxEvents=args.maxEvents)
    elif process=='sfcalc': p=PostProcessor(args.outputfile,files,cut=None, branchsel=None, outputbranchsel="keep_and_drop.txt", modules=mods,provenance=False,maxEvents=args.maxEvents)
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
    parser.add_argument('-f', '--isFastSim', action="store",  default = False,
                        help = "Input file is fastsim (Default: false)")
    parser.add_argument('-D', '--isData',  type=str, default = "",
                        help = "Data era (B, C, D, ...).  Using this flag also switches the procesor to data mode. (Default: None, i.e. MC )")
    parser.add_argument('-d', '--dataEra',    action="store",  type=str, default = "",
                        help = "Data era (B, C, D, ...).  Using this flag also switches the procesor to data mode. (Default: None, i.e. MC )")
    parser.add_argument('-s', '--sampleName',    action="store",  type=str, default = "",
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
    parser.add_argument('-p', '--process', type=str, default = '',
			help = "Type of QCD process to do (jetres or smear)")
    args = parser.parse_args()
    main(args)
