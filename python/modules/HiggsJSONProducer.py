import os, sys
from os import path
import io
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import math
import numpy as np
import gzip
import json
from array import array
from importlib import import_module

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.tools import deltaPhi, deltaR, closest
from PhysicsTools.NanoAODTools.postprocessing.framework.treeReaderArrayTools import *
from rootpy.tree import Tree, TreeModel, IntCol, FloatArrayCol

class HiggsJSONProducer(Module):
    def __init__(self, processName):
        self.metBranchName = "MET"
        self.filename=processName

    def beginJob(self):
        pass
    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.fout = gzip.open(self.filename, 'w')

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.fout.close()

    def higg2json(self, svfit, idx):
        j = [{'M':float(svfit[idx].Mass), 'pt':float(svfit[idx].Pt), 'phi':float(svfit[idx].Phi), 'eta':float(svfit[idx].Eta)}]
        j.append({'M':float(svfit[idx].tau1Mass), 'pt':float(svfit[idx].tau1Pt), 'phi':float(svfit[idx].tau1Phi), 'eta':float(svfit[idx].tau1Eta)})
        j.append({'M':float(svfit[idx].tau2Mass), 'pt':float(svfit[idx].tau2Pt), 'phi':float(svfit[idx].tau2Phi), 'eta':float(svfit[idx].tau2Eta)})
        return j

    def analyze(self, event):
        ## Getting objects
        svfit     = Collection(event, "SVFit")
        svfitmet  = Collection(event, "SVFitMET")

        for idx in xrange(len(svfit)):
            if svfit[svfit[idx].Index].Pt > 0.:
                j = self.higg2json(svfit, svfit[idx].Index)
                self.fout.write((json.dumps(j)+'\n').encode('utf-8'))

        return True
        
        
 # define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed
