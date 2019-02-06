#!/usr/bin/env python
# coding: utf-8

import os
import re
import time
import subprocess
import glob
import tarfile
import shutil
import getpass
import argparse
from collections import defaultdict

# TODO: set OutDir (and ProjectName?) to be modified based on input filelist location
DelExe    = '../Stop0l_postproc.py'
OutDir = '/store/user/%s/StopStudy' %  getpass.getuser()
tempdir = '/uscms_data/d3/%s/condor_temp/' % getpass.getuser()
ProjectName = 'PostProcess'
argument = "--inputFiles=%s.$(Process).list "
sendfiles = ["../keep_and_drop.txt"]

def tar_cmssw():
    print("Tarring up CMSSW, ignoring file larger than 100MB")
    cmsswdir = os.environ['CMSSW_BASE']
    cmsswtar = os.path.abspath('%s/CMSSW.tar.gz' % tempdir)
    if os.path.exists(cmsswtar):
        ans = raw_input('CMSSW tarball %s already exists, remove? [yn] ' % cmsswtar)
        if ans.lower()[0] == 'y':
            os.remove(cmsswtar)
        else:
            return

    def exclude(tarinfo):
        if tarinfo.size > 100*1024*1024:
            tarinfo = None
            return tarinfo
        exclude_patterns = ['/.git/', '/tmp/', '/jobs.*/', '/logs/', '/.SCRAM/', '.pyc']
        for pattern in exclude_patterns:
            if re.search(pattern, tarinfo.name):
                print('ignoring %s in the tarball', tarinfo.name)
                tarinfo = None
                break
        return tarinfo

    with tarfile.open(cmsswtar, "w:gz") as tar:
        tar.add(cmsswdir, arcname=os.path.basename(cmsswdir), filter=exclude)
    return cmsswtar

def ConfigList(config):
    process = defaultdict(dict)
    #TODO: Split between sample set and sample collection configs
    lines = open(config).readlines()
    for line in lines:
        if(len(line) <= 0 or line[0] == '#'):
            continue
        entry = line.split(",")
        stripped_entry = [ i.strip() for i in entry]
        print(stripped_entry)
        process[stripped_entry[0]] = {
            "Filepath__" : "%s/%s" % (stripped_entry[1], stripped_entry[2]),
            "isData" : len(stripped_entry) == 6,
            "isFastSim" : "fastsim" in stripped_entry[0],
            "crossSection":  float(stripped_entry[4]) * float(stripped_entry[7]),
            "nEvents":  int(stripped_entry[5])- int(stripped_entry[6]),
            "era" : 2017, #Temp
        }

    return process


# Process = {
    # "ProcessName" : ['Path to filelist', split lines]
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ SM ~~~~~
#Use SampleSet and SampleCollection configs to create a list instead.
#Then add to Process from there? How to determine the line split number? There does seem to be a default for that.
# }

def Condor_Sub(condor_file):
    curdir = os.path.abspath(os.path.curdir)
    os.chdir(os.path.dirname(condor_file))
    print "To submit condor with " + condor_file
    os.system("condor_submit " + condor_file)
    os.chdir(curdir)


def SplitPro(key, file, fraction=10):
    splitedfiles = []
    filelistdir = tempdir + '/' + "FileList"
    try:
        os.makedirs(filelistdir)
    except OSError:
        pass

    filename = os.path.abspath(file)
    if fraction == 1:
        splitedfiles.append(os.path.abspath(filename))
        shutil.copy2(os.path.abspath(filename), "%s/%s" % (filelistdir, os.path.basename(filename)))
        return splitedfiles

    print(filename)
    f = open(filename, 'r')
    lines = f.readlines()
    if len(lines) <= fraction:
        lineperfile = 1
        fraction = len(lines)
    else:
        lineperfile = len(lines) / fraction
        if len(lines) % fraction > 0:
            lineperfile += 1


    for i in range(0, fraction):
        wlines = []
        if i == fraction - 1 :
            wlines = lines[lineperfile*i :]
        else:
            wlines = lines[lineperfile*i : lineperfile*(i+1)]
        if len(wlines) > 0:
            outf = open("%s/%s.%d.list" % (filelistdir, key, i), 'w')
            outf.writelines(wlines)
            splitedfiles.append(os.path.abspath("%s/%s.%d.list" % (filelistdir, key, i)))
        outf.close()

    return splitedfiles

def my_process(args):
    ## temp dir for submit
    global tempdir
    global ProjectName
    ProjectName = time.strftime('%b%d') + ProjectName
    tempdir = tempdir + os.getlogin() + "/" + ProjectName +  "/"
    try:
        os.makedirs(tempdir)
    except OSError:
        pass

    ## Create the output directory
    outdir = OutDir +  "/" + ProjectName + "/"
    try:
        os.makedirs("/eos/uscms/%s" % outdir)
    except OSError:
        pass


    ## Update RunHT.csh with DelDir and pileups
    RunHTFile = tempdir + "/" + "RunExe.csh"
    with open(RunHTFile, "wt") as outfile:
        for line in open("RunExe.csh", "r"):
            line = line.replace("DELSCR", os.environ['SCRAM_ARCH'])
            line = line.replace("DELDIR", os.environ['CMSSW_VERSION'])
            line = line.replace("DELEXE", DelExe.split('/')[-1])
            line = line.replace("OUTDIR", outdir)
            outfile.write(line)

    ### Create Tarball
    Tarfiles = []
    NewNpro = {}

    ##Read config file
    Process = ConfigList(os.path.abspath(args.config))
    for key, sample in Process.items():
        print("Getting process: " + key + " " + sample['Filepath__'])
        npro = SplitPro(key, sample['Filepath__'])
        Tarfiles+=npro
        NewNpro[key] = len(npro)

    Tarfiles.append(os.path.abspath(DelExe))
    tarballname ="%s/%s.tar.gz" % (tempdir, ProjectName)
    with tarfile.open(tarballname, "w:gz", dereference=True) as tar:
        [tar.add(f, arcname=f.split('/')[-1]) for f in Tarfiles]
        tar.close()
    tarballname += " , %s, " % tar_cmssw()
    tarballname += " , ".join([os.path.abspath(i) for i in sendfiles])
    print(tarballname)

    ### Update condor files
    for name, sample in Process.items():
        arg = "\nArguments = --inputfile={common_name}.$(Process).list ".format(common_name=name)
        # arg = "\nArguments = --inputfile={common_name}.$(Process).list --outputfile={common_name}_$(Process).root ".format(common_name=name)
        for k, v in sample.items():
            if "__" not in k:
                arg+=" --%s=%s" % (k, v)
        arg += "\nQueue {number} \n".format(number = NewNpro[name])

        ## Prepare the condor file
        condorfile = tempdir + "/" + "condor_" + ProjectName + "_" + name
        with open(condorfile, "wt") as outfile:
            for line in open("condor_template", "r"):
                line = line.replace("EXECUTABLE", os.path.abspath(RunHTFile))
                line = line.replace("TARFILES", tarballname)
                line = line.replace("TEMPDIR", tempdir)
                line = line.replace("PROJECTNAME", ProjectName)
                #line = line.replace("OUTDIR", sample[1] + "_" + ProjectName)
                line = line.replace("ARGUMENTS", arg)
                outfile.write(line)

        Condor_Sub(condorfile)

def GetProcess(key, value):
    if len(value) == 1:
        return SplitPro(key, value[0], 1)
    else :
        return SplitPro(key, value[0], value[1])

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='NanoAOD postprocessing.')
    parser.add_argument('-c', '--config',
        default = "sampleconfig.cfg",
        help = 'Path to the input filelist.')

    args = parser.parse_args()
    my_process(args)
