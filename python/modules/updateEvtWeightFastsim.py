import os, sys
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import math

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

class UpdateEvtWeightFastsim(Module):
    def __init__(self, isData, nEvent, Process):
        self.isData = isData
        self.nEvent = nEvent
        self.process = Process
	self.xsRoot = os.environ["CMSSW_BASE"] + "/src/PhysicsTools/NanoSUSYTools/data/sigXSec/xSec.root"
	self.sTophist = "stop_xsection"
	self.gluinohist = "gluino_xsection"

    def beginJob(self):
        pass
    def endJob(self):
        pass

    def loadhisto(self,filename,hname):
        file =ROOT.TFile.Open(filename)
        hist_ = file.Get(hname)
        hist_.SetDirectory(0)
        return hist_


    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("Stop0l_evtWeight",         "F")
	self.stopXSEC=self.loadhisto(self.xsRoot,self.sTophist)
	self.gluinoXSEC=self.loadhisto(self.xsRoot,self.gluinohist)


    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
	stop0l = Object(event, "Stop0l")

        sign = 1
        if not self.isData:
            initgenWeight = getattr(event, "genWeight")
            sign          = 1 if initgenWeight > 0 else -1

	xs = 0
	if "T2" in self.process: 
		xs = self.stopXSEC.GetBinContent(self.stopXSEC.GetXaxis().FindBin(stop0l.MotherMass))
	else:
		xs = self.gluinoXSEC.GetBinContent(self.gluinoXSEC.GetXaxis().FindBin(stop0l.MotherMass))

        neweight = xs/self.nEvent * sign

        ### Store output
        self.out.fillBranch("Stop0l_evtWeight",        neweight)
        return True
