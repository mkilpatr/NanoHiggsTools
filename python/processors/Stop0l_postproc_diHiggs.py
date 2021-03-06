#!/usr/bin/env python
import os, sys
import argparse
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
from importlib import import_module
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor

from PhysicsTools.NanoHiggsTools.modules.eleMiniCutIDProducer import *
from PhysicsTools.NanoHiggsTools.modules.TauMVAObjectsProducer import *
from PhysicsTools.NanoHiggsTools.modules.updateEvtWeight import *
from PhysicsTools.NanoHiggsTools.modules.Stop0lObjectsProducer import *
from PhysicsTools.NanoHiggsTools.modules.Stop0lBaselineProducer import *
from PhysicsTools.NanoHiggsTools.modules.HiggsJSONProducer import *
from PhysicsTools.NanoHiggsTools.modules.HiggsLundVarsProducer import *
from PhysicsTools.NanoHiggsTools.modules.CutProducer import *

def main(args):
    isdata = len(args.dataEra) > 0
    isfastsim = args.isFastSim
    isVBF = args.sampleName.startswith("VBF")
    process = args.process

    mods = [eleMiniCutID(),
            UpdateEvtWeight(isdata, args.crossSection, args.nEvents, args.sampleName),
            Stop0lObjectsProducer(args.era),
            Stop0lBaselineProducer(args.era, isData=isdata, isFastSim=isfastsim),
            TauMVAObjectsProducer(isVBF=isVBF),
    ]

    if process == 'json':
        mods = [HiggsJSONProducer(args.sampleName, args.match, args.debug)]
    elif process == 'lund':
        mods = [TauMVAObjectsProducer(isVBF=isVBF), 
                HiggsLundVarsProducer(args.sampleName, args.match)]
    elif process == 'dihiggs':
        mods = [TauMVAObjectsProducer(isVBF=isVBF)]
    elif process == 'cut':
        mods = [CutProducer()]

    cut = "1 == 1"
    if args.region == "emu":      cut = "(SVFit_channel[SVFit_Index[0]] == 5 && SVFit_DZeta[SVFit_Index[0]] > -35 && SVFit_elecMuonMT[SVFit_Index[0]] < 60)"
    elif args.region == "ehad":   cut = "(SVFit_channel[SVFit_Index[0]] == 1 && SVFit_PassTight[SVFit_Index[0]] && SVFit_tau1_elecMT[SVFit_Index[0]] < 50)"
    elif args.region == "muhad":  cut = "(SVFit_channel[SVFit_Index[0]] == 0 && SVFit_PassTight[SVFit_Index[0]] && SVFit_tau1_muMT[SVFit_Index[0]] < 50)"
    elif args.region == "hadhad": cut = "(SVFit_channel[SVFit_Index[0]] == 2 && SVFit_PassTight[SVFit_Index[0]] && SVFit_ditauDR[SVFit_Index[0]] > 0.5 && SVFit_ditauPt[SVFit_Index[0]] > 50)"
    else: cut = "((SVFit_channel[SVFit_Index[0]] == 5 && SVFit_DZeta[SVFit_Index[0]] > -35 && SVFit_elecMuonMT[SVFit_Index[0]] < 60) || (SVFit_channel[SVFit_Index[0]] == 1 && SVFit_PassTight[SVFit_Index[0]] && SVFit_tau1_elecMT[SVFit_Index[0]] < 50) || (SVFit_channel[SVFit_Index[0]] == 0 && SVFit_PassTight[SVFit_Index[0]] && SVFit_tau1_muMT[SVFit_Index[0]] < 50) || (SVFit_channel[SVFit_Index[0]] == 2 && SVFit_PassTight[SVFit_Index[0]] && SVFit_ditauDR[SVFit_Index[0]] > 0.5 && SVFit_ditauPt[SVFit_Index[0]] > 50))"

    files = []
    if len(args.inputfile) > 5 and args.inputfile[0:5] == "file:":
        #This is just a single test input file
        files.append(args.inputfile[5:])
    else:
        #this is a file list
        with open(args.inputfile) as f:
            files = [line.strip() for line in f]

    if process == 'json':      p=PostProcessor(args.outputfile,files,cut="Pass_NJets30 && SVFitMET_isValid && Pass_EventFilter && Pass_JetID && nJets30 >=2 && SVFit_nPassMediumElecMuon <= 2 && " + cut, branchsel=None, outputbranchsel="keep_and_drop_train.txt", modules=mods,provenance=False,maxEvents=args.maxEvents)
    if process == 'lund':      p=PostProcessor(args.outputfile,files,cut="Pass_NJets30 && SVFitMET_isValid && Pass_EventFilter && Pass_JetID && nJets30 >=2 && SVFit_nPassMediumElecMuon <= 2", branchsel=None, outputbranchsel="keep_and_drop_train.txt", modules=mods,provenance=False,maxEvents=args.maxEvents)
    elif process == 'dihiggs': p=PostProcessor(args.outputfile,files,cut="nJet >= 2 && SVFitMET_isInteresting", branchsel=None, outputbranchsel="keep_and_drop_train.txt", modules=mods,provenance=False,maxEvents=args.maxEvents)
    elif process == 'cut':     p=PostProcessor(args.outputfile,files,cut="nJet >= 2 && SVFitMET_isInteresting", branchsel=None, outputbranchsel="keep_and_drop_train.txt", modules=mods,provenance=False,maxEvents=args.maxEvents)
    else:                      p=PostProcessor(args.outputfile,files,cut=None, branchsel=None, outputbranchsel="keep_and_drop.txt", modules=mods,provenance=False,maxEvents=args.maxEvents)
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
                        help = 'Maximum number of events to process (Default: all events)')
    parser.add_argument('-p', '--process', type=str, default = "",
                        help = "Type of QCD process to do (jetres or smear)")
    parser.add_argument('-j', '--match', type=str, default = "GenPart",
                        help = "Type of particle match for JSON files")
    parser.add_argument('-r', '--region', type=str, default = "",
                        help = "Which Region do you want to cut")
    parser.add_argument('-db', '--debug', type=bool, default = False,
                        help = "Print debug statements")
    args = parser.parse_args()
    main(args)
