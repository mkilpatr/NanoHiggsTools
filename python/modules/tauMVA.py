#!/usr/bin/env python
import os, sys
import ROOT
import math
import numpy as np
from array import array
ROOT.PyConfig.IgnoreCommandLineOptions = True
from importlib import import_module
from os import system, environ

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.tools import deltaPhi, deltaR, closest
from PhysicsTools.NanoSUSYTools.modules.xgbHelper import XGBHelper

class tauMVAProducer(Module):
    def __init__(self):
	self.writeHistFile=True
	#still trying to find appropriate cut, but the this is the best training model
	self.tauMVADisc = 0.855342
	self.bdt_file = environ["CMSSW_BASE"] + "/src/PhysicsTools/NanoSUSYTools/data/tauMVA/tauMVA-xgb_nvar13_eta0_003000_maxdepth10.model"
	self.bdt_vars = ['pt', 'abseta', 'chiso0p1', 'chiso0p2', 'chiso0p3', 'chiso0p4', 'totiso0p1', 'totiso0p2', 'totiso0p3', 'totiso0p4', 'neartrkdr', 'contjetdr', 'contjetcsv']
	self.xgb = XGBHelper(self.bdt_file, self.bdt_vars)

    def beginJob(self,histFile=None,histDirName=None):
   	pass
    def endJob(self):
	pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("taumva", 		"F", lenVar="nPFcand")
	self.out.fillBranch("TauMVA_Stop0l",    "I")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def SelTauMVA(self, mva):
	if mva > self.tauMVADisc:
		return True
	else:
		return False

    def analyze(self, event):
        ## Getting objects
	jets	  = Collection(event, "Jet")
	pfcand    = Collection(event, "PFcand")
	eventNum  = event.event

        pfchargedhads = []
	mva = {}
	mva_ = []
      
	for pfc in pfcand:
		pt 	     = min(pfc.pt,float(300.0))
		abseta       = min(abs(pfc.eta), float(2.4))
		chiso0p1     = min(pfc.chiso0p1,float(700.0))
		chiso0p2     = min(pfc.chiso0p2,float(700.0))
		chiso0p3     = min(pfc.chiso0p3,float(700.0))
		chiso0p4     = min(pfc.chiso0p4,float(700.0))
		totiso0p1    = min(pfc.totiso0p1,float(700.0))
		totiso0p2    = min(pfc.totiso0p2,float(700.0))
		totiso0p3    = min(pfc.totiso0p3,float(700.0))
		totiso0p4    = min(pfc.totiso0p4,float(700.0))
		neartrkdr    = pfc.nearestTrkDR
		jetmatch     = (pfc.contJetIndex > -1) and (jets[pfc.contJetIndex].pt >= 20.0) and (abs(jets[pfc.contJetIndex].eta) < 2.4)
		jetdr        = deltaR(jets[pfc.contJetIndex].eta, jets[pfc.contJetIndex].phi, pfc.eta, pfc.phi) if jetmatch else -1.0
		jetcsv       = jets[pfc.contJetIndex].btagDeepB if jetmatch else -1.0
		
		contjetdr  = min(float(0.4), jetdr)
		if(contjetdr < 0.0): contjetdr = 0.0
		contjetcsv =  jetcsv
		if(contjetcsv < 0.0): contjetcsv = 0.0
	
		mva = {self.bdt_vars[0]: pt, 
		       self.bdt_vars[1]: abseta,
		       self.bdt_vars[2]: chiso0p1, 
		       self.bdt_vars[3]: chiso0p2, 
		       self.bdt_vars[4]: chiso0p3, 
		       self.bdt_vars[5]: chiso0p4, 
		       self.bdt_vars[6]: totiso0p1, 
		       self.bdt_vars[7]: totiso0p2, 
		       self.bdt_vars[8]: totiso0p3, 
		       self.bdt_vars[9]: totiso0p4, 
		       self.bdt_vars[10]: neartrkdr, 
		       self.bdt_vars[11]: contjetdr, 
		       self.bdt_vars[12]: contjetcsv}
		mva_.append(self.xgb.eval(mva))
	self.TauMVA_Stop0l = map(self.SelTauMVA, mva_)

	#print "mva output: ", mva_
        self.out.fillBranch("taumva", 		mva)
	self.out.fillBranch("TauMVA_Stop0l", sum(self.TauMVA_Stop0l))
		
	return True
