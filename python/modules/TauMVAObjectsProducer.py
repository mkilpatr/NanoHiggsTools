import os, sys
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import math
import numpy as np
from array import array
from importlib import import_module

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.tools import deltaPhi, deltaR, closest
from PhysicsTools.NanoAODTools.postprocessing.framework.treeReaderArrayTools import *
from rootpy.tree import Tree, TreeModel, IntCol, FloatArrayCol

#2016 MC: https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation80XReReco#Data_MC_Scale_Factors_period_dep
#2017 MC: https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation94X

class TauMVAObjectsProducer(Module):
    def __init__(self):
        self.metBranchName = "MET"
	self.p_tauminus = 15
	self.p_Z0       = 23
	self.p_Wplus    = 24

    def beginJob(self):
        pass
    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
	self.out.branch("Pfcand_pt",		"F")
	self.out.branch("Pfcand_eta",		"F")
	self.out.branch("Pfcand_phi",		"F")
	self.out.branch("Pfcand_mass",		"F")
	self.out.branch("Pfcand_dz",		"F")
	self.out.branch("Pfcand_chiso0p1",		"F")
	self.out.branch("Pfcand_chiso0p2",		"F")
	self.out.branch("Pfcand_chiso0p3",		"F")
	self.out.branch("Pfcand_chiso0p4",		"F")
	self.out.branch("Pfcand_totiso0p1",		"F")
	self.out.branch("Pfcand_totiso0p2",		"F")
	self.out.branch("Pfcand_totiso0p3",		"F")
	self.out.branch("Pfcand_totiso0p4",		"F")
	self.out.branch("Pfcand_tackiso",		"F")
	self.out.branch("Pfcand_nearphopt",		"F")
	self.out.branch("Pfcand_nearphoeta",	"F")
	self.out.branch("Pfcand_nearphophi",	"F")
	self.out.branch("Pfcand_nearestTrkDR",	"F")
	self.out.branch("Pfcand_contJetIndex",	"I")
	self.out.branch("Pfcand_contjetdr",		"F")
	self.out.branch("Pfcand_contjetcsv",	"F")
	self.out.branch("isTau_MVA",		"O")

    def isA(self, particleID, p):
	return abs(p) == particleID

    def SelTau(self, genpart, pfc):
	for g in genpart:
		genPartMom = g.genPartIdxMother
		if self.isA(self.p_tauminus, g.pdgId) and deltaR(g.eta, g.phi, pfc.eta, pfc.phi) < 0.2 and (self.isA(self.p_Z0, genPartMom) or self.isA(self.p_Wplus, genPartMom)):
			return True
	return False

    def analyze(self, event):
        ## Getting objects
	jets	  = Collection(event, "Jet")
	genpart   = Collection(event, "GenPart")
	pfcand    = Collection(event, "PFcand")
	eventNum  = event.event

	for iP in xrange(len(pfcand)):
		pfc = pfcand[iP]
		#if(abs(pfc.eta) > 2.4 or pfc.nearestTrkDR > 2.38): continue
	
		jetmatch = (pfc.contJetIndex > -1) and (jets[pfc.contJetIndex].pt >= 20.0) and (abs(jets[pfc.contJetIndex].eta) < 2.4)
		jetdr = deltaR(jets[pfc.contJetIndex].eta, jets[pfc.contJetIndex].phi, pfc.eta, pfc.phi) if jetmatch else -1.0
		jetcsv = jets[pfc.contJetIndex].btagDeepB if jetmatch else -1.0
		self.out.fillBranch("Pfcand_pt",		pfc.pt)
		self.out.fillBranch("Pfcand_eta",		pfc.eta)
		self.out.fillBranch("Pfcand_phi",		pfc.phi)
		self.out.fillBranch("Pfcand_mass",		pfc.mass)
		self.out.fillBranch("Pfcand_dz",		pfc.dz)
		self.out.fillBranch("Pfcand_chiso0p1",	pfc.chiso0p1)
		self.out.fillBranch("Pfcand_chiso0p2",	pfc.chiso0p2)
		self.out.fillBranch("Pfcand_chiso0p3",	pfc.chiso0p3)
		self.out.fillBranch("Pfcand_chiso0p4",	pfc.chiso0p4)
		self.out.fillBranch("Pfcand_totiso0p1",	pfc.totiso0p1)
		self.out.fillBranch("Pfcand_totiso0p2",	pfc.totiso0p2)
		self.out.fillBranch("Pfcand_totiso0p3",	pfc.totiso0p3)
		self.out.fillBranch("Pfcand_totiso0p4",	pfc.totiso0p4)
		self.out.fillBranch("Pfcand_tackiso",	pfc.trackiso)
		self.out.fillBranch("Pfcand_nearphopt",	pfc.nearphopt)
		self.out.fillBranch("Pfcand_nearphoeta",	pfc.nearphoeta)
		self.out.fillBranch("Pfcand_nearphophi",	pfc.nearphophi)
		self.out.fillBranch("Pfcand_nearestTrkDR",	pfc.nearestTrkDR)
		self.out.fillBranch("Pfcand_contJetIndex",	pfc.contJetIndex)
		self.out.fillBranch("Pfcand_contjetdr",	jetdr)
		self.out.fillBranch("Pfcand_contjetcsv",	jetcsv)
		self.out.fillBranch("isTau_MVA",		self.SelTau(genpart, pfc))
		self.out.fill()

	return True


 # define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed
