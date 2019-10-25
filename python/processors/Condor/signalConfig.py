#!/usr/bin/env python
# coding: utf-8

from __future__ import division
import os
import re
import time
import subprocess
import glob
import tarfile
import shutil
import getpass
import math
import argparse
import uproot
from collections import defaultdict
from multiprocessing import Pool

#DelExe    = '../Stop0l_postproc.py'
tempdir = '/uscms_data/d3/%s/condor_temp/' % getpass.getuser()
ShortProjectName = 'PostProcess'
VersionNumber = '_v4'
argument = "--inputFiles=%s.$(Process).list "
#sendfiles = ["../keep_and_drop.txt", "../keep_and_drop_tauMVA.txt"]
sendfiles = ["../keep_and_drop.txt", "../keep_and_drop_tauMVA.txt", "../keep_and_drop_train.txt", "../keep_and_drop_LL.txt", "../keep_and_drop_res.txt", "../keep_and_drop_QCD.txt"]
TTreeName = "Events"
NProcess = 10

def ConfigList(config):
    #Allow for grabbing the era from the config file name instead. Only does so if era argument is not given.
    process = defaultdict(dict)
    #TODO: Split between sample set and sample collection configs
    conf=open("SampleSets_Postprocessed_Signal.cfg","w+")
    lines = open(config).readlines()
    for line_ in lines:
        line = line_.strip()
        if(len(line) <= 0 or line[0] == '#'):
            continue
        entry = line.split(",")
        stripped_entry = [ i.strip() for i in entry]
        #print(stripped_entry)
        replaced_outdir = stripped_entry[1].replace("Pre","Post")
        process[stripped_entry[0]] = {
            #Note that anything appended with __ will not be passed along. These are for bookkeeping. Furthermore, Outpath is not used if an output directory argument is given.
            "Filepath__" : "%s/%s" % (stripped_entry[1], stripped_entry[2]),
            #"Outpath__" : "%s" % (stripped_entry[1]) + "/" + ShortProjectName + VersionNumber + "/" + stripped_entry[0]+"/", #old
            #"Outpath__" : "%s" % (replaced_outdir) + VersionNumber + "/" + stripped_entry[0] + "/", #new
            "Outpath__" : "%s" % (replaced_outdir) + VersionNumber + "/", #new
            "isData__" : "Data" in stripped_entry[0],
            "isFastSim" : "fastsim" in stripped_entry[0], #isFastSim is a toggle in Stop0l_postproc.py, so it should be sent with no value.
        }
        if process[stripped_entry[0]]["isData__"]:
            process[stripped_entry[0]].update( {
                "dataEra": stripped_entry[0][-1], #Example naming convention: Data_MET_2018_PeriodC. Alternate option: match "Period", take location + 6.
                "crossSection":  float(stripped_entry[4]) , #storing lumi for data
                "nEvents":  int(stripped_entry[5]),
            })
        else:
            process[stripped_entry[0]].update( {
                "crossSection":  float(stripped_entry[4]) * float(stripped_entry[7]),
                "nEvents":  int(stripped_entry[5]) - int(stripped_entry[6]), # using all event weight
                "sampleName": stripped_entry[0], #process
                "totEvents__":  int(stripped_entry[5]) + int(stripped_entry[6]), # using all event weight
            })
	
	listCommand = "ls " + stripped_entry[1] + "/" + stripped_entry[0] + " > " + stripped_entry[0] + ".txt"
	print(listCommand)
	os.system(listCommand)
	
	sep = open(stripped_entry[0] + ".txt").readlines()
	for line_ in sep:
		name = os.path.splitext(line_)[0]
		f=open(name+".txt","w+")
		f.write("root://cmseos.fnal.gov/"+stripped_entry[1]+"/"+stripped_entry[0]+"/"+name+".root")
		f.close()
		conf.write(name + ", " + os.getcwd() + ", " + name+ ".txt, Events, 1.0, " + stripped_entry[5] + ", " + stripped_entry[6] + ", 1.0")
		conf.write("\n")

    conf.close()

    return process

def my_process(args):
    ##Read config file
    Process = ConfigList(os.path.abspath(args.config))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='NanoAOD postprocessing.')
    parser.add_argument('-c', '--config',
        default = "sampleconfig.cfg",
        help = 'Path to the input config file.')

    args = parser.parse_args()
    my_process(args)
