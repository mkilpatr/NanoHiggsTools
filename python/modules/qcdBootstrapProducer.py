#!/usr/bin/env python
import os, sys
import ROOT
from array import array
ROOT.PyConfig.IgnoreCommandLineOptions = True
from importlib import import_module
from os import system, environ

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.tools import deltaPhi, deltaR, closest
from PhysicsTools.NanoAODTools.postprocessing.framework.treeReaderArrayTools import *
from rootpy.tree import Tree, TreeModel, IntCol, FloatArrayCol


class qcdBootstrapProducer(Module): 
    def __init__(self):
        self.writeHistFile=True
        self.nBootstraps = 100

    def beginJob(self,histFile=None,histDirName=None):
        pass

    def endJob(self):
        pass 

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
	self.out = wrappedOutputTree
	self.out.branch("nBootstrapWeight",        "I")
	self.out.branch("bootstrapWeight",         "I", lenVar="nBootstrapWeight")

    
    def analyze(self, event):
	#Need to initialize a random seed
	ROOT.gRandom.SetSeed(123456)

	b = []        
	for iB in xrange(self.nBootstraps) :
		b.append(min(255,ROOT.gRandom.Poisson(1)))


	self.out.fillBranch("nBootstrapWeight",        self.nBootstraps)
	self.out.fillBranch("bootstrapWeight",         b)
        return True    
