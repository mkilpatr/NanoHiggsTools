#!/usr/bin/env python
import os, sys
import argparse
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
from importlib import import_module
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoSUSYTools.modules.LLObjectsProducer import *
from PhysicsTools.NanoSUSYTools.modules.TopReweightProducer import *
from PhysicsTools.NanoSUSYTools.modules.qcdSFProducer import *
from PhysicsTools.NanoSUSYTools.modules.updateEvtWeightFastsim import *
from PhysicsTools.NanoSUSYTools.modules.SoftBDeepAK8SFProducer import SoftBDeepAK8SFProducer

DataDepInputs = {
    "MC": {
        "2016" : {"bTagEff": "allInOne_bTagEff_deepCSVb_med.root",
                  "pileup_Data": "Cert271036_284044_23Sep2016ReReco_Collisions16.root",
                  "pileup_MC": "pileup_profile_2016.root",
                  "JERMC": "Summer16_25nsV1_MC",
                  "JECMC": "Summer16_07Aug2017_V11_MC",
                  "redoJEC": False,
                  "taggerWD": "TopTaggerCfg-DeepResolved_DeepCSV_GR_nanoAOD_2016_v1.0.6",
                  "nISRjets": "allInOne_ISRWeight.root",
                 },
        "2017" : {"bTagEff": "allInOne_bTagEff_deepCSVb_med.root",
                  "pileup_Data": "Cert294927_306462_EOY2017ReReco_Collisions17.root",
                  "pileup_MC": "pileup_profile_2017.root",
                  "JERMC": "Fall17_V3_MC",
                  "JECMC": "Fall17_17Nov2017_V32_MC",
                  "redoJEC": False,
                  "taggerWD": "TopTaggerCfg-DeepResolved_DeepCSV_GR_nanoAOD_2017_v1.0.6",
                  "nISRjets": "allInOne_ISRWeight.root",
                 },
        "2018" : {"bTagEff": "allInOne_bTagEff_deepCSVb_med.root",
                  "pileup_Data": "ReReco2018ABC_PromptEraD_Collisions18.root",
                  "pileup_MC": "pileup_profile_2018.root",
                  "JERMC": "Autumn18_V7_MC",
                  "JECMC": "Autumn18_V19_MC",
                  "redoJEC": True,
                  "taggerWD": "TopTaggerCfg-DeepResolved_DeepCSV_GR_nanoAOD_2018_v1.0.6",
                  "nISRjets": "allInOne_ISRWeight.root",
                 }
    },

    "FASTSIM": {
        "2016" : {"bTagEff": "FastSim2016AllSamples.root",
                  "pileup_Data": "Cert271036_284044_23Sep2016ReReco_Collisions16.root",
                  "pileup_MC": "pileup_fastsim_2016.root",
                  "JERMC": "Summer16_25nsV1_MC",
                  "JECMC": "Spring16_25nsFastSimV1_MC",
                  "redoJEC": False,
                  "taggerWD": "TopTaggerCfg-DeepResolved_DeepCSV_GR_nanoAOD_2016_v1.0.6",
                  "nISRjets": "fastsim_2016.root",
                 },
        "2017" : {"bTagEff": "FastSim2017AllSamples.root",
                  "pileup_Data": "Cert294927_306462_EOY2017ReReco_Collisions17.root",
                  "pileup_MC": "pileup_fastsim_2017.root",
                  "JERMC": "Fall17_V3_MC",
                  "JECMC": "Fall17_FastsimV1_MC",
                  "redoJEC": True,
                  "taggerWD": "TopTaggerCfg-DeepResolved_DeepCSV_GR_nanoAOD_2017_v1.0.6",
                  "nISRjets": "fastsim_2017.root",
                 },
        "2018" : {"bTagEff": "FastSim2018AllSamples.root",
                  "pileup_Data": "ReReco2018ABC_PromptEraD_Collisions18.root",
                  "pileup_MC": "pileup_fastsim_2018.root",
                  "JERMC": "Autumn18_V1_MC",
                  "JECMC": "Autumn18_FastSimV1_MC",
                  "redoJEC": True,
                  "taggerWD": "TopTaggerCfg-DeepResolved_DeepCSV_GR_nanoAOD_2018_v1.0.6",
                  "nISRjets": "fastsim_2018.root",
                 }
    },

    "Data": {
        "2016B" : { "JEC": "Summer16_07Aug2017BCD_V11_DATA",
                    "redoJEC": False,
                    "taggerWD": "TopTaggerCfg-DeepResolved_DeepCSV_GR_nanoAOD_2016_v1.0.6",
                   },
        "2016C" : { "JEC": "Summer16_07Aug2017BCD_V11_DATA",
                    "redoJEC": False,
                    "taggerWD": "TopTaggerCfg-DeepResolved_DeepCSV_GR_nanoAOD_2016_v1.0.6",
                   },
        "2016D" : { "JEC": "Summer16_07Aug2017BCD_V11_DATA",
                    "redoJEC": False,
                    "taggerWD": "TopTaggerCfg-DeepResolved_DeepCSV_GR_nanoAOD_2016_v1.0.6",
                   },
        "2016E" : { "JEC": "Summer16_07Aug2017EF_V11_DATA",
                    "redoJEC": False,
                    "taggerWD": "TopTaggerCfg-DeepResolved_DeepCSV_GR_nanoAOD_2016_v1.0.6",
                   },
        "2016F" : { "JEC": "Summer16_07Aug2017EF_V11_DATA",
                    "redoJEC": False,
                    "taggerWD": "TopTaggerCfg-DeepResolved_DeepCSV_GR_nanoAOD_2016_v1.0.6",
                   },
        "2016G" : { "JEC": "Summer16_07Aug2017GH_V11_DATA",
                    "redoJEC": False,
                    "taggerWD": "TopTaggerCfg-DeepResolved_DeepCSV_GR_nanoAOD_2016_v1.0.6",
                   },
        "2016H" : { "JEC": "Summer16_07Aug2017GH_V11_DATA",
                    "redoJEC": False,
                    "taggerWD": "TopTaggerCfg-DeepResolved_DeepCSV_GR_nanoAOD_2016_v1.0.6",
                   },

        "2017B" : { "JEC": "Fall17_17Nov2017B_V32_DATA",
                    "redoJEC": False,
                    "taggerWD": "TopTaggerCfg-DeepResolved_DeepCSV_GR_nanoAOD_2017_v1.0.6",
                   },
        "2017C" : { "JEC": "Fall17_17Nov2017C_V32_DATA",
                    "redoJEC": False,
                    "taggerWD": "TopTaggerCfg-DeepResolved_DeepCSV_GR_nanoAOD_2017_v1.0.6",
                   },
        "2017D" : { "JEC": "Fall17_17Nov2017DE_V32_DATA",
                    "redoJEC": False,
                    "taggerWD": "TopTaggerCfg-DeepResolved_DeepCSV_GR_nanoAOD_2017_v1.0.6",
                   },
        "2017E" : { "JEC": "Fall17_17Nov2017DE_V32_DATA",
                    "redoJEC": False,
                    "taggerWD": "TopTaggerCfg-DeepResolved_DeepCSV_GR_nanoAOD_2017_v1.0.6",
                   },
        "2017F" : { "JEC": "Fall17_17Nov2017F_V32_DATA",
                    "redoJEC": False,
                    "taggerWD": "TopTaggerCfg-DeepResolved_DeepCSV_GR_nanoAOD_2017_v1.0.6",
                   },

        "2018A" : { "JEC": "Autumn18_RunA_V19_DATA",
                    "redoJEC": True,
                    "taggerWD": "TopTaggerCfg-DeepResolved_DeepCSV_GR_nanoAOD_2018_v1.0.6",
                   },
        "2018B" : { "JEC": "Autumn18_RunB_V19_DATA",
                    "redoJEC": True,
                    "taggerWD": "TopTaggerCfg-DeepResolved_DeepCSV_GR_nanoAOD_2018_v1.0.6",
                   },
        "2018C" : { "JEC": "Autumn18_RunC_V19_DATA",
                    "redoJEC": True,
                    "taggerWD": "TopTaggerCfg-DeepResolved_DeepCSV_GR_nanoAOD_2018_v1.0.6",
                   },

        "2018D" : { "JEC": "Autumn18_RunD_V19_DATA",
                    "redoJEC": True,
                    "taggerWD": "TopTaggerCfg-DeepResolved_DeepCSV_GR_nanoAOD_2018_v1.0.6",
                   },
            }
}

def main(args):
    isdata = len(args.dataEra) > 0
    process = args.process
    isfastsim = args.isFastSim
    isSUSY = args.sampleName.startswith("SMS_")

    if isdata:
        dataType="Data"
        if not args.era + args.dataEra in DataDepInputs[dataType].keys():
            print "ERROR: Era \"" + args.era + "\" not recognized"
            exit(0)
    elif isfastsim:
        dataType="FASTSIM"
        if not args.era + args.dataEra in DataDepInputs[dataType].keys():
            print "ERROR: Era \"" + args.era + "\" not recognized"
            exit(0)
    else:
        dataType = "MC"
        if not args.era in DataDepInputs[dataType].keys():
            print "ERROR: Era \"" + args.era + "\" not recognized"
            exit(0)

    taggerWorkingDirectory = os.environ["CMSSW_BASE"] + "/src/PhysicsTools/NanoSUSYTools/python/processors/" + DataDepInputs[dataType][args.era if not isdata else (args.era + args.dataEra)]["taggerWD"]


    mods = [
	TopReweightProducer(args.era, args.sampleName, isData=isdata),
	#LLObjectsProducer(args.era, args.sampleName, isData=isdata),
        #SoftBDeepAK8SFProducer(args.era, taggerWorkingDirectory, isData=isdata, isFastSim=isfastsim, sampleName=args.sampleName),
    ]

#    if not isdata:
#	mods.append(LLObjectsProducer(args.era, args.sampleName, isData=isdata, applyUncert="JESUp"))
#	mods.append(LLObjectsProducer(args.era, args.sampleName, isData=isdata, applyUncert="JESDown"))
#	mods.append(LLObjectsProducer(args.era, args.sampleName, isData=isdata, applyUncert="METUnClustUp"))
#	mods.append(LLObjectsProducer(args.era, args.sampleName, isData=isdata, applyUncert="METUnClustDown"))

    if "QCD" in args.sampleName:
	mods.append(qcdSFProducer(args.era))

    if isSUSY:
	mods.append(UpdateEvtWeightFastsim(isdata, args.nEvents, args.sampleName))

    files = []
    if len(args.inputfile) > 5 and args.inputfile[0:5] == "file:":
        #This is just a single test input file
        files.append(args.inputfile[5:])
    else:
        #this is a file list
        with open(args.inputfile) as f:
            files = [line.strip() for line in f]

    if process=="limits": p=PostProcessor(args.outputfile,files,cut="Pass_MET && Pass_NJets30", branchsel=None, outputbranchsel="keep_and_drop_limits.txt", modules=mods,provenance=False,maxEvents=args.maxEvents)
    else: 		  p=PostProcessor(args.outputfile,files,cut="Pass_MET && Pass_NJets30", branchsel=None, outputbranchsel="keep_and_drop_LL.txt", modules=mods,provenance=False,maxEvents=args.maxEvents)
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
    parser.add_argument('-p', '--process', type=str, default = "",
                        help = "Type of QCD process to do (jetres or smear)")
    args = parser.parse_args()
    main(args)
