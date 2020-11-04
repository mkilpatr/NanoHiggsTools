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

    def beginJob(self):
        pass
    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
	self.out.branch("pt", "F")
	self.out.branch("mt", "F")
	self.out.branch("misset", "F")
	self.out.branch("abseta", "F")
	self.out.branch("absdz", "F")
	#self.out.branch("taumva", "F")
	self.out.branch("chiso0p1", "F")
	self.out.branch("chiso0p2", "F")
	self.out.branch("chiso0p3", "F")
	self.out.branch("chiso0p4", "F")
	self.out.branch("totiso0p1", "F")
	self.out.branch("totiso0p2", "F")
	self.out.branch("totiso0p3", "F")
	self.out.branch("totiso0p4", "F")
	self.out.branch("neartrkdr", "F")
	self.out.branch("contjetdr", "F")
	self.out.branch("contjetcsv", "F")
	self.out.branch("gentaumatch", "O")
	self.out.branch("ptmatch", "F")
	self.out.branch("etamatch", "F")
        self.out.branch("Jet_dijetMass", "F")
        self.out.branch("Tau_dijetMass", "F")
        for j in xrange(3):
            self.out.branch("Jet_pt" + str(j+1), "F")
            self.out.branch("Jet_eta" + str(j+1), "F")
            self.out.branch("Jet_phi" + str(j+1), "F")
            self.out.branch("Jet_mass" + str(j+1), "F")


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

    def addFourVec(self, obj1, obj2):
        tot = ROOT.TLorentzVector()
        v1 = ROOT.TLorentzVector()
        v2 = ROOT.TLorentzVector()
        v1.SetPtEtaPhiM(obj1.pt, obj1.eta, obj1.phi, obj1.mass)
        v2.SetPtEtaPhiM(obj2.pt, obj2.eta, obj2.phi, obj2.mass)
        tot = (v1 + v2)
        return tot

    def diJetMass(self, jet, genjet):
        j4vec = []
        for gj in xrange(len(genjet)):
            if gj == 2: break
            for j in xrange(len(jet)):
                if gj != jet[j].genJetIdx: continue
                j4vec.append(jet[j])

        if len(j4vec) == 2: return self.addFourVec(j4vec[0], j4vec[1]).M()
        elif len(j4vec) == 1: return j4vec[0].mass
        else: return -1

    def diTauMass(self, tau):
        tot = ROOT.TLorentzVector()
        tau1 = ROOT.TLorentzVector()
        tau2 = ROOT.TLorentzVector()

        if len(tau) > 1:
            tau1.SetPtEtaPhiM(tau[0].pt, tau[0].eta, tau[0].phi, tau[0].mass)
            tau2.SetPtEtaPhiM(tau[1].pt, tau[1].eta, tau[1].phi, tau[1].mass)
            tot = (tau1 + tau2)
            return tot.M()
        elif len(tau) == 1:
            return tau[0].mass
        else:
            return -1

    def analyze(self, event):
        ## Getting objects
	met	  = Object(event, "MET")
	jets	  = Collection(event, "Jet")
	genjets	  = Collection(event, "GenJet")
	genPart   = Collection(event, "GenPart")
	pfcand    = Collection(event, "Tau")
	eventNum  = event.event

	taudecayprods = []
	nGenTaus = 0
	nGenHadTaus = 0
	nGenLeptons = 0
	nGenChHads = 0
	nGenChHadsAcc = 0
	for p in genPart:
	        if p.statusFlags & 4:
			#print "staitusFlag =", p.statusFlags
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

	
	misset = met.pt
	nGenChHads = len(taudecayprods)

        for j in xrange(len(jets)):
            if j == 3: break
            self.out.fillBranch("Jet_pt"   + str(j + 1),   jets[j].pt)
            self.out.fillBranch("Jet_eta"  + str(j + 1),  jets[j].eta)
            self.out.fillBranch("Jet_phi"  + str(j + 1),  jets[j].phi)
            self.out.fillBranch("Jet_mass" + str(j + 1), jets[j].mass)

        self.out.fillBranch("Jet_dijetMass", self.diJetMass(jets, genjets))
        self.out.fillBranch("Tau_dijetMass", self.diTauMass(pfcand))

	#for pfc in pfcand:
      
	#	match = False
	#	tmpDr = 0.05
	#	kpt = 0.01
	#	ptmatch = -1.0
	#	etamatch = -10
	#	
	#	
	#	for genchhad in taudecayprods:
	#		dpt = 0.0
	#		if(genchhad.pt>0.5): 
	#			dpt = abs(1.0 - pfc.pt/genchhad.pt)
	#		if((deltaR(pfc.eta, pfc.phi, genchhad.eta, genchhad.phi) +  kpt*dpt) < tmpDr and dpt < 0.4):
	#			tmpDr = deltaR(pfc.eta, pfc.phi, genchhad.eta, genchhad.phi) + kpt*dpt
	#			match = True
	#			ptmatch = genchhad.pt
	#			etamatch = genchhad.eta
	#	
	#	if(pfc.pt > 10.0 and abs(pfc.eta) < 2.4):
	#	
	#		pt = min(pfc.pt,float(300.0))
	#		mt = self.computeMT(pfc, met, pfcand)
	#		
	#		abseta       = abs(pfc.eta)
	#		absdz        = abs(pfc.dz)
	#		chiso0p1     = min(pfc.chiso0p1,float(700.0))
	#		chiso0p2     = min(pfc.chiso0p2,float(700.0))
	#		chiso0p3     = min(pfc.chiso0p3,float(700.0))
	#		chiso0p4     = min(pfc.chiso0p4,float(700.0))
	#		totiso0p1    = min(pfc.totiso0p1,float(700.0))
	#		totiso0p2    = min(pfc.totiso0p2,float(700.0))
	#		totiso0p3    = min(pfc.totiso0p3,float(700.0))
	#		totiso0p4    = min(pfc.totiso0p4,float(700.0))
	#		neartrkdr    = pfc.nearestTrkDR
	#		jetmatch     = (pfc.contJetIndex > -1) and (jets[pfc.contJetIndex].pt >= 20.0) and (abs(jets[pfc.contJetIndex].eta) < 2.4)
	#		jetdr        = deltaR(jets[pfc.contJetIndex].eta, jets[pfc.contJetIndex].phi, pfc.eta, pfc.phi) if jetmatch else -1.0
	#		jetcsv       = jets[pfc.contJetIndex].btagDeepB if jetmatch else -1.0
	#		
	#		contjetdr  = min(float(0.4), jetdr)
	#		if(contjetdr < 0.0): contjetdr = 0.0
	#		contjetcsv =  jetcsv
	#		if(contjetcsv < 0.0): contjetcsv = 0.0
	#		if(match and nGenHadTaus > 0): gentaumatch = True
	#		else:                          gentaumatch = False
	#		
	#		self.out.fillBranch("pt",		pt)
	#		self.out.fillBranch("abseta",		abseta)
	#		self.out.fillBranch("absdz",		absdz)
	#		self.out.fillBranch("chiso0p1",		chiso0p1)
	#		self.out.fillBranch("chiso0p2",		chiso0p2)
	#		self.out.fillBranch("chiso0p3",		chiso0p3)
	#		self.out.fillBranch("chiso0p4",		chiso0p4)
	#		self.out.fillBranch("totiso0p1",	totiso0p1)
	#		self.out.fillBranch("totiso0p2",	totiso0p2)
	#		self.out.fillBranch("totiso0p3",	totiso0p3)
	#		self.out.fillBranch("totiso0p4",	totiso0p4)
	#		self.out.fillBranch("neartrkdr",	neartrkdr)
	#		self.out.fillBranch("contjetdr",	contjetdr)
	#		self.out.fillBranch("contjetcsv",	contjetcsv)
	#		self.out.fillBranch("mt", 		mt)
	#		self.out.fillBranch("misset", 		misset)
	#		self.out.fillBranch("gentaumatch", 	gentaumatch)
	#		self.out.fillBranch("ptmatch", 		ptmatch)
	#		self.out.fillBranch("etamatch", 	etamatch)
	#		self.out.fill()
		
	return True


 # define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed
