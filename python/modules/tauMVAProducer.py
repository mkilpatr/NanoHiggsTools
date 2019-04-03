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
        self.metBranchName = "MET"
	self.p_tauminus = 15
	self.p_Z0       = 23
	self.p_Wplus    = 24
	self.p_gamma    = 22
	self.pfhfhad = 1 
	self.pfem = 2 
	self.pfelectron = 11 
	self.p_nu_e = 12
	self.pfmuon = 13 
	self.p_nu_mu = 14
	self.p_nu_tau = 16
	self.pfphoton = 22 
	self.pfh0 = 130 
	self.pfhplus = 211
	self.tauMVADisc = 0.855342
	self.bdt_file_eta3 	= environ["CMSSW_BASE"] + "/src/PhysicsTools/NanoSUSYTools/data/tauMVA/tauMVA-xgb_nvar13_eta0_300000_maxdepth10.model"
	self.bdt_file_eta03 	= environ["CMSSW_BASE"] + "/src/PhysicsTools/NanoSUSYTools/data/tauMVA/tauMVA-xgb_nvar13_eta0_030000_maxdepth10.model"
	self.bdt_file_eta003 	= environ["CMSSW_BASE"] + "/src/PhysicsTools/NanoSUSYTools/data/tauMVA/tauMVA-xgb_nvar13_eta0_003000_maxdepth10.model"
	self.bdt_file_eta0003 	= environ["CMSSW_BASE"] + "/src/PhysicsTools/NanoSUSYTools/data/tauMVA/tauMVA-xgb_nvar13_eta0_000300_maxdepth10.model"
	self.bdt_file_eta00003 	= environ["CMSSW_BASE"] + "/src/PhysicsTools/NanoSUSYTools/data/tauMVA/tauMVA-xgb_nvar13_eta0_000030_maxdepth10.model"
	self.bdt_vars = ['pt', 'abseta', 'chiso0p1', 'chiso0p2', 'chiso0p3', 'chiso0p4', 'totiso0p1', 'totiso0p2', 'totiso0p3', 'totiso0p4', 'neartrkdr', 'contjetdr', 'contjetcsv']
	self.xgb_eta3 		= XGBHelper(self.bdt_file_eta3, self.bdt_vars)
	self.xgb_eta03 		= XGBHelper(self.bdt_file_eta03, self.bdt_vars)
	self.xgb_eta003 	= XGBHelper(self.bdt_file_eta003, self.bdt_vars)
	self.xgb_eta0003 	= XGBHelper(self.bdt_file_eta0003, self.bdt_vars)
	self.xgb_eta00003 	= XGBHelper(self.bdt_file_eta00003, self.bdt_vars)

    def beginJob(self,histFile=None,histDirName=None):
   	pass
    def endJob(self):
	pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
	self.out.branch("nGenHadTaus", 		"I")
	self.out.branch("nGenTaus", 		"I")
	self.out.branch("nGenChHads", 		"I")
	self.out.branch("nGenChHadsAcc", 	"I")
	self.out.branch("pt", 			"F", lenVar="nPFcand")
	self.out.branch("mt", 			"F", lenVar="nPFcand")
	self.out.branch("misset", 		"F")
	self.out.branch("abseta", 		"F", lenVar="nPFcand")
	self.out.branch("absdz", 		"F", lenVar="nPFcand")
	self.out.branch("gentaumatch", 		"O", lenVar="nPFcand")
        self.out.branch("taumva_eta3", 		"F", lenVar="nPFcand")
        self.out.branch("taumva_eta03", 	"F", lenVar="nPFcand")
        self.out.branch("taumva_eta003", 	"F", lenVar="nPFcand")
        self.out.branch("taumva_eta0003", 	"F", lenVar="nPFcand")
        self.out.branch("taumva_eta00003", 	"F", lenVar="nPFcand")
	self.out.branch("GoodTaus", 		"O", lenVar="nPFcand")
	self.out.branch("nGoodTaus", 		"I")
	self.out.branch("FakeTaus", 		"O", lenVar="nPFcand")
	self.out.branch("nFakeTaus", 		"I")
	#self.out.branch("TauMVA_Stop0l", 	"I", lenVar="nPFcand")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def SelTauMVA(self, mva):
	if mva > self.tauMVADisc:
		return True
	else:
		return False

    def isA(self, particleID, p):
	return abs(p) == particleID

    def getNearPhotonIndex(self, pfc, pfcands):
	minPhotonPt = 0.5
	maxPhotonDR = 0.2
	photonInd = -1
	maxPhotonPT = 0.0
	
	for ic in xrange(len(pfcands)):
		c = pfcands[ic]
		if(c.nearphopt < minPhotonPt): continue
		dr = deltaR(c.eta, c.phi, pfc.nearphoeta, pfc.nearphophi)
		if(dr > maxPhotonDR): continue;
		if(c.nearphopt > maxPhotonPT):
			maxPhotonPT = c.pt;
			photonInd = ic;
	
	return photonInd;

    def transverseMass(self, visible, invisible):
	cosDPhi   = np.cos( deltaPhi(visible.Phi(), invisible.phi) );
	return np.sqrt( 2 * visible.Pt() * invisible.pt * (1 - cosDPhi) );
        
    def computeMT(self, pfc, met, pfcands):
	photonInd = self.getNearPhotonIndex(pfc, pfcands);
	candP4 = ROOT.TLorentzVector()
	candP4.SetPtEtaPhiM(pfc.pt, pfc.eta, pfc.phi, pfc.mass)
	if(photonInd > -1): 
		pfcand_buff = ROOT.TLorentzVector()
		pfcand_buff.SetPtEtaPhiM(pfcands[photonInd].pt, pfcands[photonInd].eta, pfcands[photonInd].phi, pfcands[photonInd].mass)
		candP4+=pfcand_buff;
	return self.transverseMass(candP4, met);

    def analyze(self, event):
        ## Getting objects
	met	  = Object(event, self.metBranchName)
	jets	  = Collection(event, "Jet")
	genPart   = Collection(event, "GenVisTau")
	pfcand    = Collection(event, "PFcand")
	eventNum  = event.event

        pfchargedhads = []
        pfphotons = []
	mva = {}
	mva_eta3 = []
	mva_eta03 = []
	mva_eta003 = []
	mva_eta0003 = []
	mva_eta00003 = []
	GoodTaus_ = []
	FakeTaus_ = []
      
        #for c in pfcand :
	#	if c.pdgId == self.pfhplus:  pfchargedhads.append(c)
	#	if c.pdgId == self.pfphoton: pfphotons.append(c)
      
        taudecayprods = [];
	nGenHadTaus = 0
	nGenTaus = 0
	nGenChHadsAcc = 0
        for p in genPart :
		if(self.isA(self.p_Z0, p.genPartIdxMother) or self.isA(self.p_Wplus, p.genPartIdxMother)):
			nGenTaus += 1
			if p.status != 15:
				taudecayprods.append(p)
				nGenChHadsAcc += 1
			nGenHadTaus += 1
      
        nGenChHads = len(taudecayprods)
      
        misset = met.pt

	mt_ = []
	gentaumatch_ = [] 
	for pfc in pfcand:
      
		match = False
		tmpDr = 0.05
		kpt = 0.01
		ptmatch = -1.0
		etamatch = -10
		GoodTaus = False
		FakeTaus = False
		
		for genchhad in taudecayprods:
			dpt = 0.0
			if(genchhad.pt>0.5): dpt = abs(1.0 - pfc.pt/genchhad.pt);
			if((deltaR(pfc.eta, pfc.phi, genchhad.eta, genchhad.phi) +  kpt*dpt) < tmpDr and dpt < 0.4):
			  tmpDr = deltaR(pfc.eta, pfc.phi, genchhad.eta, genchhad.phi) + kpt*dpt
			  match = True
			  ptmatch = genchhad.pt
			  etamatch = genchhad.eta
		
		if(pfc.pt > 10.0 and abs(pfc.eta) < 2.4):
		
			pt = min(pfc.pt,float(300.0))
			mt = self.computeMT(pfc, met, pfcand)
			mt_.append(mt)
			
			abseta       = abs(pfc.eta)
			absdz        = abs(pfc.dz)
			#taumva       = pfc.taudisc;
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
			if(match and nGenHadTaus > 0): gentaumatch = True
			else:                          gentaumatch = False
			gentaumatch_.append(gentaumatch)		
	
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
			mva_eta3.append(self.xgb_eta3.eval(mva))
			mva_eta03.append(self.xgb_eta03.eval(mva))
			mva_eta003.append(self.xgb_eta003.eval(mva))
			mva_eta0003.append(self.xgb_eta0003.eval(mva))
			mva_eta00003.append(self.xgb_eta00003.eval(mva))
			if gentaumatch==1 and nGenHadTaus>0  and len(jets)>3 and misset>150 and mt<100 and pt>10 and ptmatch > 6. and absdz<0.2: GoodTaus = True
			if gentaumatch==0 and nGenHadTaus==0 and len(jets)>3 and misset>150 and mt<100 and pt>10 and absdz<0.2: FakeTaus = True
			GoodTaus_.append(GoodTaus)
			FakeTaus_.append(FakeTaus)
	#self.TauMVA_Stop0l = map(self.SelTauMVA, mva_)

	#print "mva output: ", mva_
	self.out.fillBranch("nGenHadTaus", 	nGenHadTaus)
	self.out.fillBranch("nGenTaus", 	nGenTaus)
	self.out.fillBranch("nGenChHads", 	nGenChHads)
	self.out.fillBranch("nGenChHadsAcc", 	nGenChHadsAcc)
	self.out.fillBranch("mt", 		mt_)
	self.out.fillBranch("misset", 		misset)
	self.out.fillBranch("gentaumatch", 	gentaumatch_)
        self.out.fillBranch("taumva_eta3", 	mva_eta3)
        self.out.fillBranch("taumva_eta03", 	mva_eta03)
        self.out.fillBranch("taumva_eta003", 	mva_eta003)
        self.out.fillBranch("taumva_eta0003", 	mva_eta0003)
        self.out.fillBranch("taumva_eta00003", 	mva_eta00003)
	self.out.fillBranch("GoodTaus", GoodTaus_)
	self.out.fillBranch("nGoodTaus", sum(GoodTaus_))
	self.out.fillBranch("FakeTaus", FakeTaus_)
	self.out.fillBranch("nFakeTaus", sum(FakeTaus_))
	#self.out.fillBranch("TauMVA_Stop0l", sum(self.TauMVA_Stop0l))
		
	return True
