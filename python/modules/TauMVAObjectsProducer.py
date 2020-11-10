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
	self.pfhneut = 111

    def beginJob(self):
        pass
    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
	self.out.branch("nGenHadTaus", 	 "I")
	self.out.branch("nGenTaus", 	 "I")
	self.out.branch("nGenChHads", 	 "I")
	self.out.branch("nGenLeptons", 	 "I")
        self.out.branch("nJets20",       "I")
        self.out.branch("Jet_dijetMass", "F")
        self.out.branch("Tau_dijetMass", "F")
        self.out.branch("Jet_deltaR",    "F")
        self.out.branch("Tau_deltaR",    "F")
        self.out.branch("HiggsCand_pt",  "F") 
        self.out.branch("HiggsCand_eta", "F") 
        self.out.branch("HiggsCand_phi", "F") 
        self.out.branch("HiggsCand_mass","F") 
        for j in xrange(3):
            self.out.branch("Jet_pt" + str(j+1), "F")
            self.out.branch("Jet_eta" + str(j+1), "F")
            self.out.branch("Jet_phi" + str(j+1), "F")
            self.out.branch("Jet_mass" + str(j+1), "F")


    # HAS BIT
    def hasBit(self, value,bit):
        """Check if i'th bit is set to 1, i.e. binary of 2^(i-1),
        from the right to the left, starting from position i=0."""
        # https://cms-nanoaod-integration.web.cern.ch/integration/master-102X/mc102X_doc.html#GenPart
        # Gen status flags, stored bitwise, are:
        #    0: isPrompt,                          8: fromHardProcess,
        #    1: isDecayedLeptonHadron,             9: isHardProcessTauDecayProduct,
        #    2: isTauDecayProduct,                10: isDirectHardProcessTauDecayProduct,
        #    3: isPromptTauDecayProduct,          11: fromHardProcessBeforeFSR,
        #    4: isDirectTauDecayProduct,          12: isFirstCopy,
        #    5: isDirectPromptTauDecayProduct,    13: isLastCopy,
        #    6: isDirectHadronDecayProduct,       14: isLastCopyBeforeFSR
        #    7: isHardProcess,
        ###return bin(value)[-bit-1]=='1'
        ###return format(value,'b').zfill(bit+1)[-bit-1]=='1'
        return (value & (1 << bit))>0

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
        if jet.pt < 20 or math.fabs(jet.eta) > 2.4:
            return False
        return True

    def DeltaR(self, obj):
        if len(obj) > 1: dr = deltaR(obj[0].eta, obj[0].phi, obj[1].eta, obj[1].phi)
        elif len(obj) == 1: dr = deltaR(obj[0].eta, obj[0].phi, 0, 0)
        else: dr = -1
        return dr

    def analyze(self, event):
        ## Getting objects
	met	  = Object(event, self.metBranchName)
	jets	  = Collection(event, "Jet")
	genPart   = Collection(event, "GenPart")
	pfcand    = Collection(event, "PFCands")
	eventNum  = event.event

        taudecayprods = []
        antitaudecayprods = []
	nGenTaus = 0
	nGenHadTaus = 0
	nGenLeptons = 0
	nGenChHads = 0
        decayType = [0, 0]
	for p in genPart:
	    if self.hasBit(p.statusFlags, 3):
	        nGenTaus+=1
	        lepdecay = False
                print("Event Num: {0}, Mother pdgID {1}, daughter {2}".format(eventNum, genPart[p.genPartIdxMother].pdgId, p.pdgId))
                genMomPdgId = genPart[p.genPartIdxMother].pdgId
	        if (self.isA(self.pfelectron, p.pdgId) or self.isA(self.pfmuon, p.pdgId)):
	            lepdecay = True
	            if genMomPdgId == self.p_tauminus: 
                        taudecayprods.append(p)
                        decayType[0] = lepdecay
	            elif abs(genMomPdgId) == self.p_tauminus: 
                        antitaudecayprods.append(p)
                        decayType[1] = lepdecay
	        if (not self.isA(self.p_nu_e, p.pdgId)) and (not self.isA(self.p_nu_mu, p.pdgId)):
	            if (self.isA(self.pfhplus, p.pdgId) or self.isA(321, p.pdgId) or self.isA(self.pfhneut, p.pdgId)):
	                if genMomPdgId == self.p_tauminus: 
                            taudecayprods.append(p)
                            decayType[0] = lepdecay
	                elif abs(genMomPdgId) == self.p_tauminus: 
                            antitaudecayprods.append(p)
                            decayType[1] = lepdecay
	        if not lepdecay:
	            nGenHadTaus+=1
	        if self.isA(self.pfelectron, p.pdgId) or self.isA(self.pfmuon, p.pdgId):
	            nGenLeptons+=1

        vec = [ROOT.TLorentzVector(), ROOT.TLorentzVector()]
        hvec = ROOT.TLorentzVector()
        taudr = 0.
        print("nPFCand {0}".format(len(pfcand)))
        for p in pfcand:
            trk = ROOT.TLorentzVector()

	    match = False
	    tmpDr = 0.05
	    kpt = 0.01
	    ptmatch = -1.0
	    etamatch = -10
	    
	    for genchhad in taudecayprods:
	        dpt = 0.0
	        if(genchhad.pt>0.5): 
	            dpt = abs(1.0 - p.pt/genchhad.pt)
                print("Tau Size {2} pdgId: {3} tau DeltaR: {0} + {1}".format(deltaR(p.eta, p.phi, genchhad.eta, genchhad.phi), kpt*dpt, len(taudecayprods), genchhad.pdgId))
	        if((deltaR(p.eta, p.phi, genchhad.eta, genchhad.phi) +  kpt*dpt) < tmpDr):
	            tmpDr = deltaR(p.eta, p.phi, genchhad.eta, genchhad.phi) + kpt*dpt
	            match = True
	            ptmatch = genchhad.pt
	            etamatch = genchhad.eta
                    trk.SetPtEtaPhiM(p.pt, p.eta, p.phi, p.mass)
            vec[0] += (trk)
            trk.SetPtEtaPhiM(0, 0, 0, 0)

	    for genchhad in antitaudecayprods:
	        dpt = 0.0
	        if(genchhad.pt>0.5): 
	            dpt = abs(1.0 - p.pt/genchhad.pt)
                print("AntiTau Size {2} pdgId: {3} tau DeltaR: {0} + {1}".format(deltaR(p.eta, p.phi, genchhad.eta, genchhad.phi), kpt*dpt, len(antitaudecayprods), genchhad.pdgId))
	        if((deltaR(p.eta, p.phi, genchhad.eta, genchhad.phi) +  kpt*dpt) < tmpDr):
	            tmpDr = deltaR(p.eta, p.phi, genchhad.eta, genchhad.phi) + kpt*dpt
	            match = True
	            ptmatch = genchhad.pt
	            etamatch = genchhad.eta
                    trk.SetPtEtaPhiM(p.pt, p.eta, p.phi, p.mass)
            vec[1] += (trk)
            print("pdgId {8} vec[0]: ({0}, {1}, {2}, {3}), vec[1]: ({4}, {5}, {6}, {7})".format(vec[0].Pt(), vec[0].Eta(), vec[0].Phi(), vec[0].M(), vec[1].Pt(), vec[1].Eta(), vec[1].Phi(), vec[1].M(), p.pdgId))

            #if p.pt < 20 or abs(p.eta) > 2.4: continue

        taudr = deltaR(vec[0].Eta(), vec[0].Phi(), vec[1].Eta(), vec[1].Phi())
        taumass = (vec[0] + vec[1])
        self.out.fillBranch("Tau_dijetMass", taumass.M())
        self.out.fillBranch("Tau_deltaR", taudr)

        #Add met to di-tau 4-vector
        hvec += (vec[0] + vec[1])
        metvec = ROOT.TLorentzVector()
        metvec.SetPtEtaPhiM(met.pt, 0, met.phi, 0)
        hvec += (metvec)

        self.out.fillBranch("HiggsCand_pt",   hvec.Pt())
        self.out.fillBranch("HiggsCand_eta",  hvec.Eta())
        self.out.fillBranch("HiggsCand_phi",  hvec.Phi())
        self.out.fillBranch("HiggsCand_mass", hvec.M())

        self.out.fillBranch("nGenHadTaus", 	nGenHadTaus)
	self.out.fillBranch("nGenTaus", 	nGenTaus)
	self.out.fillBranch("nGenChHads", 	nGenChHads)
	self.out.fillBranch("nGenLeptons", 	nGenLeptons)

        #Jet Variables
        self.Jet_Stop0l      = map(self.SelJets, jets)
        self.out.fillBranch("nJets20",          sum(self.Jet_Stop0l))
        jetvec = ROOT.TLorentzVector()
        jet = []
        for j in xrange(len(jets)):
            if j == 3: break
            if self.Jet_Stop0l[j]:
                jet.append(jets[j])
                self.out.fillBranch("Jet_pt"   + str(j + 1),   jets[j].pt)
                self.out.fillBranch("Jet_eta"  + str(j + 1),  jets[j].eta)
                self.out.fillBranch("Jet_phi"  + str(j + 1),  jets[j].phi)
                self.out.fillBranch("Jet_mass" + str(j + 1), jets[j].mass)

        self.out.fillBranch("Jet_dijetMass", self.addFourVec(jet).M())
        self.out.fillBranch("Jet_deltaR", self.DeltaR(jet))

		
	return True


 # define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed
