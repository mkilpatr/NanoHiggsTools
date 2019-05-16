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
from PhysicsTools.NanoSUSYTools.modules.PDFUncertaintyProducer import *
from PhysicsTools.NanoSUSYTools.modules.GenPartFilter import GenPartFilter
from PhysicsTools.NanoAODTools.postprocessing.modules.common.puWeightProducer import *
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jecUncertainties import jecUncertProducer
from PhysicsTools.NanoSUSYTools.modules.tauMVAProducer import *
from PhysicsTools.NanoSUSYTools.modules.TauMVAObjectsProducer import *
from PhysicsTools.NanoSUSYTools.modules.LLObjectsProducer import *

# JEC files are those recomended here (as of Mar 1, 2019)
# https://twiki.cern.ch/twiki/bin/view/CMS/JECDataMC#Recommended_for_MC
# Actual text files are found here
# JEC: https://github.com/cms-jet/JECDatabase/tree/master/textFiles
# JER: https://github.com/cms-jet/JRDatabase/tree/master/textFiles
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
    isdata = len(args.dataEra) > 0
    isfastsim = args.isFastSim
    process = args.process
    isfakemva = True
    iseff = True if process == "taumvacompare" else False

    if isdata and isfastsim:
        print "ERROR: It is impossible to have a dataset that is both data and fastsim"
        exit(0)

    mods = []
    if process == "train":
	mods.append(TauMVAObjectsProducer())
    elif process == "taumva" or process == "taumvacompare":
    	mods += [
    	    eleMiniCutID(),
    	    Stop0lObjectsProducer(args.era),
    	    DeepTopProducer(args.era),
    	    Stop0lBaselineProducer(args.era, isData=isdata, isFastSim=isfastsim),
	    UpdateEvtWeight(isdata, args.crossSection, args.nEvents, args.sampleName),
    	    tauMVAProducer(isFakeMVA=isfakemva, isEff=iseff, isData=isdata),
	    LLObjectsProducer(args.era),
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
    	        #puWeightProducer("auto", pufile, "pu_mc","pileup", verbose=False),
    	        # statusFlag 0x2100 corresponds to "isLastCopy and fromHardProcess"
    	        # statusFlag 0x2080 corresponds to "IsLastCopy and isHardProcess"
    	        GenPartFilter(statusFlags = [0x2100, 0x2080]),
    	    ]

    files = ["root://cmseos.fnal.gov//eos/uscms/store/user/lpcsusyhad/Stop_production/Autumn18_102X_v1/PreProcessed_22March2019/MET//2018_Data_Run2018A-17Sep2018-v1/190330_215429/0000/prod2018DATA_NANO_1-41.root"]
    #files = []
    #if len(args.inputfile) > 5 and args.inputfile[0:5] == "file:":
    #    #This is just a single test input file
    #    files.append(args.inputfile[5:])
    #else:
    #    #this is a file list
    #    with open(args.inputfile) as f:
    #        files = [line.strip() for line in f]

    if process=="train":    
	p=PostProcessor(args.outputfile,files,cut="Pass_MET & Pass_Baseline", branchsel=None, outputbranchsel="keep_and_drop_train.txt", typeofprocess="tau", modules=mods,provenance=False)
    elif process=="taumva": 
	p=PostProcessor(args.outputfile,files,cut="MET_pt > 150 & nJet > 3", branchsel=None, outputbranchsel="keep_and_drop_tauMVA.txt", modules=mods,provenance=False)
    elif process == "taumvacompare":
	p=PostProcessor(args.outputfile,files,cut="MET_pt > 150", branchsel=None, outputbranchsel="keep_and_drop_tauMVA.txt", modules=mods,provenance=False)
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
                        help = 'Maximum number of events to process (Default: all events)')
    parser.add_argument('-p', '--process', type=str, default = "",
                        help = "Type of QCD process to do (jetres or smear)")
    args = parser.parse_args()
    main(args)
