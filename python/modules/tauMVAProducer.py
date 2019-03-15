#!/usr/bin/env python
import os, sys
import ROOT
import math
ROOT.PyConfig.IgnoreCommandLineOptions = True
from importlib import import_module

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.tools import deltaPhi, deltaR, closest
from PhysicsTools.NanoSUSYTools.modules.xgbHelper import XGBHelper

class tauMVAProducer(Module):
    def __init__(self):
	self.writeHistFile=True
	self.p_tauminus = 15
	self.p_Z0       = 23
	self.p_Wplus    = 24
	self.tauMVADisc = 0.855342
	#self.bdt_file = "/eos/uscms/store/user/mkilpatr/13TeV/tauMVA/xgboost.xml"
	self.bdt_file = "/uscms_data/d3/mkilpatr/CMSSW_10_2_9/src/TauMVATraining/tauMVA-xgb.model"
	self.bdt_vars = ["Pfcand_pt", "Pfcand_eta", "Pfcand_chiso0p1", "Pfcand_chiso0p2", "Pfcand_chiso0p3", "Pfcand_chiso0p4", "Pfcand_totiso0p1", "Pfcand_totiso0p2", "Pfcand_totiso0p3", "Pfcand_totiso0p4", "Pfcand_nearestTrkDR", "Pfcand_contjetdr", "Pfcand_contjetcsv"]
	self.xgb = XGBHelper(self.bdt_file, self.bdt_vars)

    def beginJob(self,histFile=None,histDirName=None):
   	pass
    def endJob(self):
	pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("TauMVA_disc", "F", lenVar="nPFcand")
	self.out.branch("TauMVA_match", "O", lenVar="nPFcand")
	self.out.branch("TauMVA_Stop0l", "I");

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def SelTauMVA(self, mva):
	if mva > self.tauMVADisc:
		return True
	else:
		return False

    def isA(self, particleID, p):
        return abs(p) == particleID

    def SelTau(self, genpart, pfc):
        for g in genpart:
                genPartMom = g.genPartIdxMother
                if self.isA(self.p_tauminus, g.pdgId) and deltaR(g.eta, g.phi, pfc.eta, pfc.phi) < 0.2 and (self.isA(self.p_Z0, genPartMom) or self.isA(self.p_Wplus, genPartMom)):
                        return True
        return False

    def analyze(self, event):
	pfcand    = Collection(event, "PFcand")
	jet	  = Collection(event, "Jet")
	genpart   = Collection(event, "GenPart")

	mva = {}
	mva_ = []
	matchTau_ = []
	for pfc in pfcand :
		pfcIdx = pfc.contJetIndex
		jetmatch = pfcIdx > -1 and jet[pfcIdx].pt > 30.0 and math.fabs(jet[pfcIdx].eta) < 2.4;
		contjetdr = deltaR(pfc, jet[pfcIdx]) if jetmatch else -1.0;
		contjetcsv = jet[pfcIdx].btagDeepB if jetmatch else -1.0;
		contjetdr = min(float(0.4), contjetdr)
		mva = {self.bdt_vars[0]: min(pfc.pt, float(300.0)), 
		       self.bdt_vars[1]: min(abs(pfc.eta), float(2.4)), 
		       self.bdt_vars[2]: min(pfc.chiso0p1 ,float(700)), 
		       self.bdt_vars[3]: min(pfc.chiso0p2 ,float(700)), 
		       self.bdt_vars[4]: min(pfc.chiso0p3 ,float(700)), 
		       self.bdt_vars[5]: min(pfc.chiso0p4 ,float(700)), 
		       self.bdt_vars[6]: min(pfc.totiso0p1,float(700)), 
		       self.bdt_vars[7]: min(pfc.totiso0p2,float(700)), 
		       self.bdt_vars[8]: min(pfc.totiso0p3,float(700)), 
		       self.bdt_vars[9]: min(pfc.totiso0p4,float(700)), 
		       self.bdt_vars[10]: min(pfc.nearestTrkDR, float(2.38)), 
		       self.bdt_vars[11]: 0.0 if contjetdr < 0.0 else contjetdr, 
		       self.bdt_vars[12]: 0.0 if contjetdr < 0.0 else contjetcsv}
		mva_.append(self.xgb.eval(mva))
		matchTau_.append(self.SelTau(genpart, pfc))

	self.TauMVA_Stop0l = map(self.SelTauMVA, mva_)
	#print "mva output: ", mva_
	self.out.fillBranch("TauMVA_disc", mva_)
	self.out.fillBranch("TauMVA_match", matchTau_)
	self.out.fillBranch("TauMVA_Stop0l", sum(self.TauMVA_Stop0l))

        return True
