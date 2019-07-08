#!/usr/bin/env python
import os, sys
import ROOT
import math
import numpy as np
import xgboost as xgb
import logging
from array import array
ROOT.PyConfig.IgnoreCommandLineOptions = True
from importlib import import_module
from os import system, environ

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.tools import deltaPhi, deltaR, closest

class XGBHelper:
    def __init__(self, model_file, var_list):
        self.bst = xgb.Booster(params={'nthread': 1}, model_file=model_file)
        self.var_list = var_list
        logging.info('Load XGBoost model %s, input variables:\n  %s' % (model_file, str(var_list)))

    def eval(self, inputs):
        dmat = xgb.DMatrix(np.array([[inputs[k] for k in self.var_list]]), feature_names=self.var_list)
        return self.bst.predict(dmat)[0]

class tauMVAProducer(Module):
    def __init__(self, isFakeMVA = False, isEff = False, isData = False):
	self.writeHistFile=True
	self.isFakeMVA = isFakeMVA 
	self.isEff = isEff
	self.isData = isData
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
	self.pfhneu = 111
	self.bdt_file_eta3 	= environ["CMSSW_BASE"] + "/src/PhysicsTools/NanoSUSYTools/data/tauMVA/tauMVA-xgb_nvar13_eta0_300000_maxdepth10.model"
	self.bdt_file_eta03 	= environ["CMSSW_BASE"] + "/src/PhysicsTools/NanoSUSYTools/data/tauMVA/tauMVA-xgb_nvar13_eta0_030000_maxdepth10.model"
	self.bdt_file_eta003 	= environ["CMSSW_BASE"] + "/src/PhysicsTools/NanoSUSYTools/data/tauMVA/tauMVA-xgb_nvar13_eta0_003000_maxdepth10.model"
	self.bdt_file_eta0003 	= environ["CMSSW_BASE"] + "/src/PhysicsTools/NanoSUSYTools/data/tauMVA/tauMVA-xgb_nvar13_eta0_000300_maxdepth10.model"
	self.bdt_file_eta00003 	= environ["CMSSW_BASE"] + "/src/PhysicsTools/NanoSUSYTools/data/tauMVA/tauMVA-xgb_nvar13_eta0_000030_maxdepth10.model"
	self.bdt_file 		= environ["CMSSW_BASE"] + "/src/PhysicsTools/NanoSUSYTools/data/tauMVA/tauMVA-xgb_nvar13_eta0_030000_maxdepth10.model"
	self.bdt_vars = ['pt', 'abseta', 'chiso0p1', 'chiso0p2', 'chiso0p3', 'chiso0p4', 'totiso0p1', 'totiso0p2', 'totiso0p3', 'totiso0p4', 'neartrkdr', 'contjetdr', 'contjetcsv']
	self.xgb 		= XGBHelper(self.bdt_file, self.bdt_vars)
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
	self.out.branch("nGenEle", 		"I")
	self.out.branch("nGenMuon", 		"I")
	self.out.branch("nGenLeptons", 		"I")
	self.out.branch("pt", 			"F", lenVar="PFcand")
	self.out.branch("mt", 			"F", lenVar="PFcand")
	self.out.branch("misset", 		"F")
	self.out.branch("abseta", 		"F", lenVar="PFcand")
	self.out.branch("absdz", 		"F", lenVar="PFcand")
	self.out.branch("gentaumatch", 		"O", lenVar="PFcand")
	self.out.branch("GoodTaus", 		"O", lenVar="PFcand")
	self.out.branch("nGoodTaus", 		"I")
	self.out.branch("FakeTaus", 		"O", lenVar="PFcand")
	self.out.branch("nFakeTaus", 		"I")
	self.out.branch("TauCR", 		"O", lenVar="PFcand")
	self.out.branch("nTauCR", 		"I")
	if self.isEff:
		self.out.branch("taumva",               	"F", lenVar="nPFcand")
		self.out.branch("TauMVA_Stop0l_71",     	"O", lenVar="nPFcand")
		self.out.branch("nTauMVA_71",           	"I")
		self.out.branch("TauMVA_Stop0l_73",     	"O", lenVar="nPFcand")
		self.out.branch("nTauMVA_73",           	"I")
		self.out.branch("TauMVA_Stop0l_71_pt10to20",    "O", lenVar="nPFcand")
		self.out.branch("nTauMVA_71_pt10to20",          "I")
		self.out.branch("TauMVA_Stop0l_73_pt10to20",    "O", lenVar="nPFcand")
		self.out.branch("nTauMVA_73_pt10to20",          "I")
		self.out.branch("TauMVA_Stop0l_71_ptgeq20",     "O", lenVar="nPFcand")
		self.out.branch("nTauMVA_71_ptgeq20",           "I")
		self.out.branch("TauMVA_Stop0l_73_ptgeq20",     "O", lenVar="nPFcand")
		self.out.branch("nTauMVA_73_ptgeq20",           "I")
	else:
        	self.out.branch("taumva_eta3", 		"F", lenVar="PFcand")
        	self.out.branch("taumva_eta03", 	"F", lenVar="PFcand")
        	self.out.branch("taumva_eta003", 	"F", lenVar="PFcand")
        	self.out.branch("taumva_eta0003", 	"F", lenVar="PFcand")
        	self.out.branch("taumva_eta00003", 	"F", lenVar="PFcand")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def SelTauMVA(self, mva, tauMVADisc):
        if mva < tauMVADisc:
                return False
        return True

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
		if(dr > maxPhotonDR): continue
		if(c.nearphopt > maxPhotonPT):
			maxPhotonPT = c.pt
			photonInd = ic
	
	return photonInd

    def transverseMass(self, visible, invisible):
	cosDPhi   = np.cos( deltaPhi(visible.Phi(), invisible.phi) )
	return np.sqrt( 2 * visible.Pt() * invisible.pt * (1 - cosDPhi) )
        
    def computeMT(self, pfc, met, pfcands):
	photonInd = self.getNearPhotonIndex(pfc, pfcands)
	candP4 = ROOT.TLorentzVector()
	candP4.SetPtEtaPhiM(pfc.pt, pfc.eta, pfc.phi, pfc.mass)
	if(photonInd > -1): 
		pfcand_buff = ROOT.TLorentzVector()
		pfcand_buff.SetPtEtaPhiM(pfcands[photonInd].pt, pfcands[photonInd].eta, pfcands[photonInd].phi, pfcands[photonInd].mass)
		candP4+=pfcand_buff
	return self.transverseMass(candP4, met)

    def analyze(self, event):
        ## Getting objects
	met	  = Object(event, self.metBranchName)
	jets	  = Collection(event, "Jet")
	if not self.isData:
		genPart   = Collection(event, "GenPart")
		genVisTau = Collection(event, "GenVisTau")
	pfcand    = Collection(event, "PFcand")
	eventNum  = event.event

	taudecayprods = []
	nGenTaus = 0
	nGenHadTaus = 0
	nGenLeptons = 0
	nGenEle = 0
	nGenMuon = 0
	nGenChHads = 0
	nGenChHadsAcc = 0
	if not self.isData:
		for p in genPart:
		        if p.statusFlags & 4:
				nGenTaus+=1
		                lepdecay = False
		                if self.isA(self.pfelectron, p.pdgId) or self.isA(self.pfmuon, p.pdgId):
		                        lepdecay = True
		                        continue
		                if (not self.isA(self.p_nu_e, p.pdgId)) and (not self.isA(self.p_nu_mu, p.pdgId)):
					if (self.isA(self.pfhplus, p.pdgId) or self.isA(321, p.pdgId)):
		                                taudecayprods.append(p)
		                                if p.pt > 10.0 and abs(p.eta) < 2.4: nGenChHadsAcc+=1
		                if not lepdecay:
					nGenHadTaus+=1
				if self.isA(self.pfelectron, p.pdgId) or self.isA(self.pfmuon, p.pdgId):
					nGenLeptons+=1
			if self.isA(self.pfelectron, p.pdgId) and p.pt > 5 and abs(p.eta) < 2.5:
				nGenEle+=1
			if self.isA(self.pfmuon, p.pdgId) and p.pt > 5 and abs(p.eta) < 2.4:
				nGenMuon+=1
	

	misset = met.pt
	nGenChHads = len(taudecayprods)

        pfchargedhads = []
        pfphotons = []
	mva = {}
	mva_ = []
	mva_eta3_ = []
	mva_eta03_ = []
	mva_eta003_ = []
	mva_eta0003_ = []
	mva_eta00003_ = []
	GoodTaus_ = []
	FakeTaus_ = []
	TauCR_ = []
	mt_ = []
	pt_ = []
	abseta_ = []
	absdz_ = []
	gentaumatch_ = [] 
	for pfc in pfcand:
      
		match = False
		tmpDr = 0.05
		kpt = 0.01
		ptmatch = -1.0
		etamatch = -10
		GoodTaus = False
		FakeTaus = False
		TauCR = False
		gentaumatch = False
		mva_eta3 = -10.0
		mva_eta03 = -10.0
		mva_eta003 = -10.0
		mva_eta0003 = -10.0
		mva_eta00003 = -10.0
		mva_buff = -10.0
		mt = 999
		pt = 0.0
		abseta = 10.0
		absdz  = 10.0
		
		for genchhad in taudecayprods:
			dpt = 0.0
			if(genchhad.pt>0.5): 
				dpt = abs(1.0 - pfc.pt/genchhad.pt)
			if((deltaR(pfc.eta, pfc.phi, genchhad.eta, genchhad.phi) +  kpt*dpt) < tmpDr and dpt < 0.4):
				tmpDr = deltaR(pfc.eta, pfc.phi, genchhad.eta, genchhad.phi) + kpt*dpt
				match = True
				ptmatch = genchhad.pt
				etamatch = genchhad.eta
		
		if((pfc.pt > 10.0 and abs(pfc.eta) < 2.4 and abs(pfc.dz) < 0.2)):
			mt = self.computeMT(pfc, met, pfcand)
			if mt < 100:
				pt = min(pfc.pt,float(300.0))
				abseta       = abs(pfc.eta)
				absdz        = abs(pfc.dz)
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
				gentaumatch_.append(gentaumatch)		

				baseline = len(jets)>3 and misset>150 and mt<100 and pt>10 and absdz<0.2	
				if (self.isFakeMVA == False or self.isEff) and gentaumatch==True and nGenLeptons==0 and nGenTaus==nGenHadTaus and nGenHadTaus > 0 and baseline and  ptmatch > 6.: 
					GoodTaus = True
				if (self.isFakeMVA == True or self.isEff) and gentaumatch==False and nGenLeptons==0 and nGenTaus==0 and baseline:
					FakeTaus = True
				if self.isEff and baseline:
					TauCR = True

				if FakeTaus == True or GoodTaus == True or TauCR == True:
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
					if TauCR == True:
						mva_buff = self.xgb.eval(mva)
					else:
						mva_eta3	 = self.xgb_eta3.eval(mva)
						mva_eta03	 = self.xgb_eta03.eval(mva)
						mva_eta003	 = self.xgb_eta003.eval(mva)
						mva_eta0003	 = self.xgb_eta0003.eval(mva)
						mva_eta00003	 = self.xgb_eta00003.eval(mva)

		#print "fastsim: %d, FakeTaus: %d, GoodTaus: %d" %(self.isFakeMVA, FakeTaus, GoodTaus)
		mt_.append(mt)
		pt_.append(pt)
		abseta_.append(abseta)
		absdz_.append(absdz)
		mva_.append(mva_buff)
		mva_eta3_.append(mva_eta3)
		mva_eta03_.append(mva_eta03)
		mva_eta003_.append(mva_eta003)
		mva_eta0003_.append(mva_eta0003)
		mva_eta00003_.append(mva_eta00003)
		GoodTaus_.append(GoodTaus)
		FakeTaus_.append(FakeTaus)
		TauCR_.append(TauCR)

	cut_71 = []
	cut_73 = []
	for i in mva_:
		cut_71.append(0.71)
		cut_73.append(0.73)

        TauMVA_Stop0l_71 = map(self.SelTauMVA, mva_, cut_71)
        TauMVA_Stop0l_73 = map(self.SelTauMVA, mva_, cut_73)
	self.TauMVA_Stop0l_71_pt10to20 = [t and pt > 10 and pt < 20  for t, pt in zip(TauMVA_Stop0l_71, pt_)]
	self.TauMVA_Stop0l_73_pt10to20 = [t and pt > 10 and pt < 20  for t, pt in zip(TauMVA_Stop0l_73, pt_)]
	self.TauMVA_Stop0l_71_ptgeq20 = [t and pt > 20 for t, pt in zip(TauMVA_Stop0l_71, pt_)]
	self.TauMVA_Stop0l_73_ptgeq20 = [t and pt > 20 for t, pt in zip(TauMVA_Stop0l_73, pt_)]

	#print "mva output: ", mva_
	self.out.fillBranch("nGenHadTaus", 	nGenHadTaus)
	self.out.fillBranch("nGenTaus", 	nGenTaus)
	self.out.fillBranch("nGenChHads", 	nGenChHads)
	self.out.fillBranch("nGenChHadsAcc", 	nGenChHadsAcc)
	self.out.fillBranch("nGenEle", 		nGenEle)
	self.out.fillBranch("nGenMuon", 	nGenMuon)
	self.out.fillBranch("nGenLeptons", 	nGenLeptons)
	self.out.fillBranch("pt",               pt_)
	self.out.fillBranch("mt", 		mt_)
	self.out.fillBranch("abseta", 		abseta_)
	self.out.fillBranch("absdz", 		absdz_)
	self.out.fillBranch("misset", 		misset)
	self.out.fillBranch("gentaumatch", 	gentaumatch_)
	self.out.fillBranch("GoodTaus", GoodTaus_)
	self.out.fillBranch("nGoodTaus", sum(GoodTaus_))
	self.out.fillBranch("FakeTaus", FakeTaus_)
	self.out.fillBranch("nFakeTaus", sum(FakeTaus_))
	self.out.fillBranch("TauCR", TauCR_)
	self.out.fillBranch("nTauCR", sum(TauCR_))
	if self.isEff or self.isData:
		self.out.fillBranch("taumva",           mva_)
		self.out.fillBranch("TauMVA_Stop0l_71", TauMVA_Stop0l_71)
		self.out.fillBranch("nTauMVA_71",       sum(TauMVA_Stop0l_71))
		self.out.fillBranch("TauMVA_Stop0l_73", TauMVA_Stop0l_73)
		self.out.fillBranch("nTauMVA_73",       sum(TauMVA_Stop0l_73))
		self.out.fillBranch("TauMVA_Stop0l_71_pt10to20", self.TauMVA_Stop0l_71_pt10to20)
		self.out.fillBranch("nTauMVA_71_pt10to20",       sum(self.TauMVA_Stop0l_71_pt10to20))
		self.out.fillBranch("TauMVA_Stop0l_73_pt10to20", self.TauMVA_Stop0l_73_pt10to20)
		self.out.fillBranch("nTauMVA_73_pt10to20",       sum(self.TauMVA_Stop0l_73_pt10to20))
		self.out.fillBranch("TauMVA_Stop0l_71_ptgeq20", self.TauMVA_Stop0l_71_ptgeq20)
		self.out.fillBranch("nTauMVA_71_ptgeq20",       sum(self.TauMVA_Stop0l_71_ptgeq20))
		self.out.fillBranch("TauMVA_Stop0l_73_ptgeq20", self.TauMVA_Stop0l_73_ptgeq20)
		self.out.fillBranch("nTauMVA_73_ptgeq20",       sum(self.TauMVA_Stop0l_73_ptgeq20))
        else:
		self.out.fillBranch("taumva_eta3", 	mva_eta3_)
        	self.out.fillBranch("taumva_eta03", 	mva_eta03_)
        	self.out.fillBranch("taumva_eta003", 	mva_eta003_)
        	self.out.fillBranch("taumva_eta0003", 	mva_eta0003_)
        	self.out.fillBranch("taumva_eta00003", 	mva_eta00003_)
		
	return True
