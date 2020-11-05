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
        self.out.branch("Jet_dijetMass", "F")
        self.out.branch("Tau_dijetMass", "F")
        self.out.branch("HiggsCand_pt",  "F") 
        self.out.branch("HiggsCand_eta", "F") 
        self.out.branch("HiggsCand_phi", "F") 
        self.out.branch("HiggsCand_mass","F") 
        for j in xrange(3):
            self.out.branch("Jet_pt" + str(j+1), "F")
            self.out.branch("Jet_eta" + str(j+1), "F")
            self.out.branch("Jet_phi" + str(j+1), "F")
            self.out.branch("Jet_mass" + str(j+1), "F")


    def isA(self, particleID, p):
	return abs(p) == particleID

    def addFourVec(self, obj):
        tot = ROOT.TLorentzVector()
        v1 = ROOT.TLorentzVector()
        v2 = ROOT.TLorentzVector()
        if(len(obj) > 0): v1.SetPtEtaPhiM(obj[0].pt, obj[0].eta, obj[0].phi, obj[0].mass)
        else: v1.SetPtEtaPhiM(0, 0, 0, 0)
        if(len(obj) > 1): v2.SetPtEtaPhiM(obj[1].pt, obj[1].eta, obj[1].phi, obj[1].mass)
        else: v2.SetPtEtaPhiM(0, 0, 0, 0)
        tot = (v1 + v2)
        return tot

    def SelJets(self, jet):
        if jet.pt < 30 or math.fabs(jet.eta) > 2.4:
            return False
        return True

    def SelTauPOG(self, tau, met):
        if tau.pt < 20 or abs(tau.eta) > 2.4 or not self.isA(tau.pdgId, self.p_tauminus):
                return False
        return True

    def SelIsotrack(self, isk, met):
        iso = isk.pfRelIso03_chg
        if abs(isk.pdgId) == 11 or abs(isk.pdgId) == 13:
            if isk.pt < 5 or iso > 0.2:
                return False
        if abs(isk.pdgId) == 211:
            if isk.pt < 10 or iso > 0.1:
                return False
        return True

    def diObjMass(self, obj, isGood):
        vec = []

        for p in xrange(len(obj)):
            if isGood[p]:
                vec.append(obj[p])

        ##print "staitusFlag =", p.statusFlags
	#    nGenTaus+=1
	#    lepdecay = False
	#    if self.isA(self.pfelectron, p.pdgId) or self.isA(self.pfmuon, p.pdgId):
	#        lepdecay = True
	#        continue
	#    if (not self.isA(self.p_nu_e, p.pdgId)) and (not self.isA(self.p_nu_mu, p.pdgId)):
	#        if (self.isA(self.pfhplus, p.pdgId) or self.isA(321, p.pdgId)):
	#            taudecayprods.append(p)
	#            if p.pt > 10.0 and abs(p.eta) < 2.4: nGenChHadsAcc+=1
	#    if not lepdecay:
	#        nGenHadTaus+=1
	#    if self.isA(self.pfelectron, p.pdgId) or self.isA(self.pfmuon, p.pdgId):
	#        nGenLeptons+=1


        return self.addFourVec(vec).M()

    def analyze(self, event):
        ## Getting objects
	met	  = Object(event, self.metBranchName)
	jets	  = Collection(event, "Jet")
	genjets	  = Collection(event, "GenJet")
	genPart   = Collection(event, "GenPart")
	pfcand    = Collection(event, "PFCands")
        isotrack  = Collection(event, "IsoTrack")
	eventNum  = event.event

        self.Tau_Stop0l      = map(lambda x : self.SelTauPOG(x, met), pfcand)
        self.Jet_Stop0l      = map(self.SelJets, jets)
        self.IsoTrack_Stop0l = map(lambda x : self.SelIsotrack(x, met), isotrack)

        hvec = ROOT.TLorentzVector()
        trk = ROOT.TLorentzVector()

        for iso in xrange(len(isotrack)):
            p = isotrack[iso]
            if not self.IsoTrack_Stop0l[iso]: continue
            
            if self.isA(self.pfelectron, p.pdgId) or self.isA(self.pfmuon, p.pdgId):
                trk.SetPtEtaPhiM(p.pt, p.eta, p.phi, 0)
            elif self.isA(self.pfhplus, p.pdgId):
                trk.SetPtEtaPhiM(p.pt, p.eta, p.phi, 0)
            hvec += (trk)

        metvec = ROOT.TLorentzVector()
        metvec.SetPtEtaPhiM(met.pt, 0, met.phi, 0)
        hvec += (metvec)

        for j in xrange(len(jets)):
            if j == 3: break
            if self.Jet_Stop0l[j]:
                self.out.fillBranch("Jet_pt"   + str(j + 1),   jets[j].pt)
                self.out.fillBranch("Jet_eta"  + str(j + 1),  jets[j].eta)
                self.out.fillBranch("Jet_phi"  + str(j + 1),  jets[j].phi)
                self.out.fillBranch("Jet_mass" + str(j + 1), jets[j].mass)

        self.out.fillBranch("HiggsCand_pt",   hvec.Pt())
        self.out.fillBranch("HiggsCand_eta",  hvec.Eta())
        self.out.fillBranch("HiggsCand_phi",  hvec.Phi())
        self.out.fillBranch("HiggsCand_mass", hvec.M())
        self.out.fillBranch("Jet_dijetMass", self.diObjMass(jets,   self.Jet_Stop0l))
        self.out.fillBranch("Tau_dijetMass", self.diObjMass(pfcand, self.Tau_Stop0l))

		
	return True


 # define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed
