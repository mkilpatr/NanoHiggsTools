import os, sys
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import math
import numpy as np
from array import array
from importlib import import_module

from PhysicsTools.NanoSUSYTools.modules.HelperFunctions import *
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.tools import deltaPhi, deltaR, closest
from PhysicsTools.NanoAODTools.postprocessing.framework.treeReaderArrayTools import *
from rootpy.tree import Tree, TreeModel, IntCol, FloatArrayCol

LeptonPts = {
   0 : (23, 30),
   1 : (26, 30),
   2 : (40, 40),
   3 : (10, 10),
   4 : (13, 13),
   5 : (13, 10)
}

LeptonEtas = {
   0 : (2.1, 2.3),
   1 : (2.1, 2.3),
   2 : (2.1, 2.1),
   3 : (2.4, 2.4),
   4 : (2.5, 2.5),
   5 : (2.5, 2.4)
}

TauIds = {
   'Tight' : 5,
   'Medium' : 4,
   'Loose' : 3,
   'VLoose' : 2,
   'VVLoose' : 1
}

ElecIds = {
   'Tight' : 2,
   'Medium' : 1,
   'Loose' : 0,
}

MuonIds = {
   'HighPt' : 4,
   'Tight' : 3,
   'Medium' : 2,
   'Soft' : 1,
   'Loose' : 0,
}

pairType = {
   'kMuHad' : 0,
   'kEHad' : 1,
   'kHadHad' : 2,
   'kMuMu' : 3,    #shouldn't be needed
   'kEE' : 4,      #shouldn't be needed
   'kEMu' : 5,
   'kEEPrompt' : 6,#shouldn't be needed
   'kMuMuPrompt' : 7,#shouldn't be needed
   'kOther' : 8    # for e.g. h->bb
}

class TauMVAObjectsProducer(Module):
    def __init__(self, isVBF):
        self.metBranchName = "MET"
        self.isVBF = isVBF
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
        self.out.branch("nJets30",       "I")
        self.out.branch("SVFit_Index",       "I", lenVar="nSVFit")
        self.out.branch("SVFit_PassBaseline", "O", lenVar="nSVFit")
        self.out.branch("SVFit_PassLepton", "O", lenVar="nSVFit")
        self.out.branch("SVFit_nPassTight", "I")
        self.out.branch("SVFit_nPassMedium", "I")
        self.out.branch("SVFit_nPassVLoose", "I")
        self.out.branch("SVFit_nPassVVLoose", "I")
        self.out.branch("SVFit_nPassTightElecMuon",   "I")
        self.out.branch("SVFit_nPassMediumElecMuon",  "I")
        self.out.branch("SVFit_nPassLooseElecMuon",   "I")
        self.out.branch("SVFit_PassTight", "O", lenVar="nSVFit")
        self.out.branch("SVFit_PassMedium", "O", lenVar="nSVFit")
        self.out.branch("SVFit_PassVLoose", "O", lenVar="nSVFit")
        self.out.branch("SVFit_PassVVLoose", "O", lenVar="nSVFit")
        self.out.branch("SVFit_PassTightElecMuon",   "O", lenVar="nSVFit")
        self.out.branch("SVFit_PassMediumElecMuon",  "O", lenVar="nSVFit")
        self.out.branch("SVFit_PassLooseElecMuon",   "O", lenVar="nSVFit")
        self.out.branch("SVFit_Boosted", "O", lenVar="nSVFit")
        self.out.branch("SVFit_VBF", "O")
        self.out.branch("SVFit_elecMuonMT",  "F", lenVar="nSVFit")
        self.out.branch("SVFit_MaxMT",       "F", lenVar="nSVFit")
        self.out.branch("SVFit_tau1_elecMT", "F", lenVar="nSVFit")
        self.out.branch("SVFit_tau1_muMT",   "F", lenVar="nSVFit")
        self.out.branch("SVFit_tau1_hadMT",  "F", lenVar="nSVFit")
        self.out.branch("SVFit_tau2_muMT",   "F", lenVar="nSVFit")
        self.out.branch("SVFit_tau2_hadMT",  "F", lenVar="nSVFit")
        self.out.branch("SVFit_ditauMass", "F", lenVar="nSVFit")
        self.out.branch("SVFit_ditauPt",   "F", lenVar="nSVFit")
        self.out.branch("SVFit_ditauDR",    "F", lenVar="nSVFit")
        self.out.branch("SVFit_deltaREMu", "F", lenVar="nSVFit")
        self.out.branch("SVFit_bj1Pt",   "F")
        self.out.branch("SVFit_bj1Eta",   "F")
        self.out.branch("SVFit_bj1Phi",   "F")
        self.out.branch("SVFit_bj1Mass",   "F")
        self.out.branch("SVFit_bj1DeepCSV",   "F")
        self.out.branch("SVFit_bj2Pt",   "F")
        self.out.branch("SVFit_bj2Eta",   "F")
        self.out.branch("SVFit_bj2Phi",   "F")
        self.out.branch("SVFit_bj2Mass",   "F")
        self.out.branch("SVFit_bj2DeepCSV",   "F")
        self.out.branch("SVFit_j1Pt",   "F")
        self.out.branch("SVFit_j1Eta",   "F")
        self.out.branch("SVFit_j1Phi",   "F")
        self.out.branch("SVFit_j1Mass",   "F")
        self.out.branch("SVFit_j1DeepCSV",   "F")
        self.out.branch("SVFit_j2Pt",   "F")
        self.out.branch("SVFit_j2Eta",   "F")
        self.out.branch("SVFit_j2Phi",   "F")
        self.out.branch("SVFit_j2Mass",   "F")
        self.out.branch("SVFit_j2DeepCSV",   "F")
        self.out.branch("SVFit_j3Pt",   "F")
        self.out.branch("SVFit_j3Eta",   "F")
        self.out.branch("SVFit_j3Phi",   "F")
        self.out.branch("SVFit_j3Mass",   "F")
        self.out.branch("SVFit_j3DeepCSV",   "F")
        self.out.branch("SVFit_2tau2jetPt", "F", lenVar="nSVFit")
        self.out.branch("SVFit_dijetPt",   "F")
        self.out.branch("SVFit_dijetMass", "F")
        self.out.branch("SVFit_dibjetDR",  "F")
        self.out.branch("SVFit_MaxbjetTauDR", "F", lenVar="nSVFit")
        self.out.branch("SVFit_MinbjetTauDR", "F", lenVar="nSVFit")
        self.out.branch("SVFit_dijetDR",    "F")
        self.out.branch("SVFit_dijetDEta",  "F")
        self.out.branch("SVFit_HT",        "F")
        self.out.branch("SVFit_DZeta",        "F", lenVar="nSVFit")


    def SelJets(self, jet):
        if jet.pt < 30 or math.fabs(jet.eta) > 4.7:
            return False
        return True

    def CalMtWTau1(self, obj, met):
        if obj.tau1Pt > -100:
            #print("met: {0} tau1pt: {1} metPhi: {2} tau1Phi: {3}, dPhi: {4}, Cos: {5}".format(met.pt, obj.tau1Pt, met.phi, obj.tau1Phi, ROOT.TVector2.Phi_mpi_pi(met.phi-obj.tau1Phi), math.cos(ROOT.TVector2.Phi_mpi_pi(met.phi-obj.tau1Phi))))
            return math.sqrt( 2 * met.pt * obj.tau1Pt * (1 - math.cos(ROOT.TVector2.Phi_mpi_pi(met.phi-obj.tau1Phi))))
        else: return obj.tau1Pt

    def CalMtWTau2(self, obj, met):
        if obj.tau2Pt > -100:
            #print("met: {0} tau2pt: {1} metPhi: {2} tau2Phi: {3}, dPhi: {4}, Cos: {5}".format(met.pt, obj.tau2Pt, met.phi, obj.tau2Phi, ROOT.TVector2.Phi_mpi_pi(met.phi-obj.tau2Phi), math.cos(ROOT.TVector2.Phi_mpi_pi(met.phi-obj.tau2Phi))))
            return math.sqrt( 2 * met.pt * obj.tau2Pt * (1 - math.cos(ROOT.TVector2.Phi_mpi_pi(met.phi-obj.tau2Phi))))
        else: return obj.tau2Pt

    def CalMtW(self, obj, met):
        return math.sqrt( 2 * met.pt * obj.Pt() * (1 - math.cos(ROOT.TVector2.Phi_mpi_pi(met.phi-obj.Phi()))))

    def CalHT(self, jets):
        HT = sum([j.pt for i, j in enumerate(jets) if self.Jets_Stop0l[i]])
        return HT

    def DeltaR(self, obj):
        if len(obj) > 1: dr = deltaR(obj[0].eta, obj[0].phi, obj[1].eta, obj[1].phi)
        elif len(obj) == 1: dr = deltaR(obj[0].eta, obj[0].phi, 0, 0)
        else: dr = -1
        return dr

    def DeltaEta(self, obj):
        if len(obj) > 1: dr = obj[0].eta - obj[1].eta
        elif len(obj) == 1: dr = obj[0].eta
        else: dr = -99
        return dr

    def isChannel(self, obj, kType):
        return pairType[kType] == obj.channel


    def SelTaus(self, obj):
        if obj.channel < -100: return False
        if math.fabs(obj.tau1Eta) > LeptonEtas[obj.channel][0] or obj.tau1Pt < LeptonPts[obj.channel][0] or math.fabs(obj.tau2Eta) > LeptonEtas[obj.channel][1] or obj.tau2Pt < LeptonPts[obj.channel][1]:
            return False
        return True

    def SelLeptons(self, obj, met):
        if self.isChannel(obj, 'kMuHad') or self.isChannel(obj, 'kEHad'):
            if (obj.tau1Pt + obj.tau2Pt + met.pt) < 50:
                return False

        return True

    def SelTausID(self, obj, ID = "Medium"):
        if self.isChannel(obj, 'kEMu'): return True
        if self.isChannel(obj, 'kMuHad') or self.isChannel(obj, 'kEHad'):
            if not hasBit(obj.tau2IDjet, TauIds[ID]):
                return False
        elif self.isChannel(obj, 'kHadHad'):
            if not hasBit(obj.tau1IDjet, TauIds[ID]) or not hasBit(obj.tau2IDjet, TauIds[ID]):
                return False
        return True

    def SelElecMuonID(self, obj, ID = "Medium"):
        #if obj.channel >= 0: print("channel {0}, id1: {1}, id2: {2}".format(obj.channel, id1, id2))
        if self.isChannel(obj, 'kMuHad') and (obj.tau1Isojet >= 0.15 or not hasBit(obj.tau1IDjet, MuonIds[ID])): 
            return False
        elif self.isChannel(obj, 'kEHad') and (obj.tau1Isojet >= 0.15 or not hasBit(obj.tau1IDjet, ElecIds[ID])): 
            return False
        elif self.isChannel(obj, 'kEMu') and (obj.tau1Isojet >= 0.15 or obj.tau2Isojet >= 0.2 or not hasBit(obj.tau1IDjet, ElecIds[ID]) or not hasBit(obj.tau2IDjet, MuonIds[ID])):
            return False
        return True

    def DZeta(self, obj, met):
        if not self.isChannel(obj, 'kEMu'): return 0.
        lep1 = ROOT.TLorentzVector()
        lep2 = ROOT.TLorentzVector()
        lep1.SetPtEtaPhiM(1, 0, obj.tau1Phi, 0)
        lep2.SetPtEtaPhiM(1, 0, obj.tau2Phi, 0)
        bisector = (lep1 + lep2).Phi()
        pVis = (obj.tau1Pt + obj.tau2Pt)*np.cos(bisector)
        pMiss = met.pt*np.cos(bisector)

        return pMiss - 1.85*pVis


    def minMaxDR(self, obj1, obj2):
        dr = []
        for j2 in obj2:
            for j1 in obj1:
                dr.append(deltaR(j2.Eta(), j2.Phi(), j1.eta, j1.phi))

            if len(obj1) == 0:
                dr.append(deltaR(j2.Eta(), j2.Phi(), 0., 0.))

        return dr

    def GetSVFitSortedIdx(self, svfit):
        ptlist = []
        tightlist = []
        mediumlist = []
        vlooselist = []
        vvlooselist = []
        for sv in xrange(len(svfit)):
            ptlist.append(-svfit[sv].Pt)
            tightlist.append(-self.SVFit_Stop0l_Tight[sv])
            mediumlist.append(-self.SVFit_Stop0l_Medium[sv])
            vlooselist.append(-self.SVFit_Stop0l_VLoose[sv])
            vvlooselist.append(-self.SVFit_Stop0l_VVLoose[sv])
        sortIdx = np.lexsort((tightlist, ptlist))
        return sortIdx

    def analyze(self, event):
        ## Getting objects
        met	  = Object(event, self.metBranchName)
        jets	  = Collection(event, "Jet")
        pfcand    = Collection(event, "PFCands")
        svfit     = Collection(event, "SVFit")
        svfitmet  = Collection(event, "SVFitMET")
        genpart   = Collection(event, "GenPart")

        #Jet Variables
        self.Jets_Stop0l      = map(self.SelJets, jets)
        self.out.fillBranch("nJets30",          sum(self.Jets_Stop0l))
        HT = self.CalHT(jets)
        self.out.fillBranch("SVFit_HT",       HT)
        jetvec = ROOT.TLorentzVector()
        jet1 = ROOT.TLorentzVector()
        jet2 = ROOT.TLorentzVector()
        jet = []
        bjet = []
        for j in xrange(len(jets)):
            if self.Jets_Stop0l[j]: jet.append(jets[j])
            if jets[j].btagStop0l: bjet.append(jets[j])

        self.SVIndex = []        
        self.SVFit_Stop0l = []
        self.SVFit_Stop0l_leptonID = []
        self.SVFit_Stop0l_Tight = []
        self.SVFit_Stop0l_Medium = []
        self.SVFit_Stop0l_VLoose = []
        self.SVFit_Stop0l_VVLoose = []
        self.SVFit_Stop0l_Tight_ElecMuon = []
        self.SVFit_Stop0l_Medium_ElecMuon = []
        self.SVFit_Stop0l_Loose_ElecMuon = []
        self.SVFit_tau1_MtW = []
        self.SVFit_tau2_MtW = []
        self.SVFit_DZeta = []
        self.boosted = []
        self.drEMu = []
        self.taudr = []
        self.MT_elecMu   = [] 
        self.MT_tau1_ele = []
        self.MT_tau1_mu  = []
        self.MT_tau1_tau = []
        self.MT_tau2_mu  = []
        self.MT_tau2_tau = []
        self.MT_max      = []
        self.tausMass = []
        self.ditauPt = []
        self.MaxdrBJetTaus = []
        self.MindrBJetTaus = []
        self.SVFit_2tau2jetPt = []
        for sv in xrange(len(svfit)):
            self.SVFit_Stop0l.append( self.SelTaus(svfit[sv]))
            self.SVFit_Stop0l_leptonID.append( self.SelLeptons(svfit[sv], met))
            self.SVFit_Stop0l_Tight.append( self.SelTausID(svfit[sv], "Tight"))
            self.SVFit_Stop0l_Medium.append( self.SelTausID(svfit[sv], "Medium"))
            self.SVFit_Stop0l_VLoose.append( self.SelTausID(svfit[sv], "VLoose"))
            self.SVFit_Stop0l_VVLoose.append( self.SelTausID(svfit[sv], "VVLoose"))
            self.SVFit_Stop0l_Tight_ElecMuon.append(self.SelElecMuonID(svfit[sv], "Tight"))
            self.SVFit_Stop0l_Medium_ElecMuon.append(self.SelElecMuonID(svfit[sv], "Medium"))
            self.SVFit_Stop0l_Loose_ElecMuon.append(self.SelElecMuonID(svfit[sv], "Loose"))
            self.SVFit_tau1_MtW.append( self.CalMtWTau1(svfit[sv], met))
            self.SVFit_tau2_MtW.append( self.CalMtWTau2(svfit[sv], met))
            self.SVFit_DZeta.append( self.DZeta(svfit[sv], met))
            self.boosted.append(svfit[sv].Pt > 400.)

            self.drEMu.append(deltaR(svfit[sv].tau1Eta, svfit[sv].tau1Phi, svfit[sv].tau2Eta, svfit[sv].tau2Phi) if self.isChannel(svfit[sv], 'kEMu') else -999)
            taus = ROOT.TLorentzVector()
            tausVec = []
            tau1 = ROOT.TLorentzVector()
            tau2 = ROOT.TLorentzVector()
            MET  = ROOT.TLorentzVector()
            tau1.SetPtEtaPhiM(svfit[sv].tau1Pt, svfit[sv].tau1Eta, svfit[sv].tau1Phi, svfit[sv].tau1Mass)
            tau2.SetPtEtaPhiM(svfit[sv].tau2Pt, svfit[sv].tau2Eta, svfit[sv].tau2Phi, svfit[sv].tau2Mass)
            taus = (tau1 + tau2)
            tausVec.append(tau1)
            tausVec.append(tau2)
            self.taudr.append(deltaR(svfit[sv].tau1Eta, svfit[sv].tau1Phi, svfit[sv].tau2Eta, svfit[sv].tau2Phi))
            self.MT_elecMu.append(self.CalMtW(taus, met) if self.isChannel(svfit[sv], 'kEMu') else -999)
            self.MT_tau1_ele.append(self.SVFit_tau1_MtW[sv] if self.isChannel(svfit[sv], 'kEMu') or self.isChannel(svfit[sv], 'kEHad') else -999)
            self.MT_tau1_mu.append(self.SVFit_tau1_MtW[sv] if self.isChannel(svfit[sv], 'kMuHad') else -999)
            self.MT_tau1_tau.append(self.SVFit_tau1_MtW[sv] if self.isChannel(svfit[sv], 'kHadHad') else -999)
            self.MT_tau2_mu.append(self.SVFit_tau2_MtW[sv] if self.isChannel(svfit[sv], 'kEMu') else -999)
            self.MT_tau2_tau.append(self.SVFit_tau2_MtW[sv] if (self.isChannel(svfit[sv], 'kMuHad') or self.isChannel(svfit[sv], 'kEHad') or self.isChannel(svfit[sv], 'kHadHad')) else -999)
            self.tausMass.append(taus.M())
            MET.SetPtEtaPhiM(met.pt, 0., met.phi, 0.)
            taus = (tau1 + tau2 + MET)
            self.ditauPt.append(taus.Pt())
            drBJetTaus = self.minMaxDR(bjet, tausVec)
            drBJetTaus.sort()
            self.MaxdrBJetTaus.append(drBJetTaus[-1])
            self.MindrBJetTaus.append(drBJetTaus[0])
            jetvec = (jet1 + jet2 + tau1 + tau2 + MET)
            self.SVFit_2tau2jetPt.append(jetvec.Pt())

        self.SVIndex = self.GetSVFitSortedIdx(svfit)

        self.out.fillBranch("SVFit_Index",          self.SVIndex)
        self.out.fillBranch("SVFit_PassLepton",     self.SVFit_Stop0l_leptonID)
        self.out.fillBranch("SVFit_PassBaseline",   self.SVFit_Stop0l)
        self.out.fillBranch("SVFit_nPassTight",     sum(self.SVFit_Stop0l_Tight))
        self.out.fillBranch("SVFit_nPassMedium",    sum(self.SVFit_Stop0l_Medium))
        self.out.fillBranch("SVFit_nPassVLoose",    sum(self.SVFit_Stop0l_VLoose))
        self.out.fillBranch("SVFit_nPassVVLoose",   sum(self.SVFit_Stop0l_VVLoose))
        self.out.fillBranch("SVFit_nPassTightElecMuon",   sum(self.SVFit_Stop0l_Tight_ElecMuon))
        self.out.fillBranch("SVFit_nPassMediumElecMuon",   sum(self.SVFit_Stop0l_Medium_ElecMuon))
        self.out.fillBranch("SVFit_nPassLooseElecMuon",   sum(self.SVFit_Stop0l_Loose_ElecMuon))
        self.out.fillBranch("SVFit_PassTight",     self.SVFit_Stop0l_Tight)
        self.out.fillBranch("SVFit_PassMedium",     self.SVFit_Stop0l_Medium)
        self.out.fillBranch("SVFit_PassVLoose",     self.SVFit_Stop0l_VLoose)
        self.out.fillBranch("SVFit_PassVVLoose",    self.SVFit_Stop0l_VVLoose)
        self.out.fillBranch("SVFit_PassTightElecMuon",   self.SVFit_Stop0l_Tight_ElecMuon)
        self.out.fillBranch("SVFit_PassMediumElecMuon",  self.SVFit_Stop0l_Medium_ElecMuon)
        self.out.fillBranch("SVFit_PassLooseElecMuon",   self.SVFit_Stop0l_Loose_ElecMuon)
        self.out.fillBranch("SVFit_Boosted",   self.boosted)
        self.out.fillBranch("SVFit_VBF",   self.isVBF)
        self.out.fillBranch("SVFit_DZeta", self.SVFit_DZeta)

        self.MT_max = [max(t1, t2) for t1, t2 in zip(self.MT_tau1_ele, self.MT_tau2_mu)]
        self.out.fillBranch("SVFit_elecMuonMT", self.MT_elecMu)
        self.out.fillBranch("SVFit_MaxMT", self.MT_max)
        self.out.fillBranch("SVFit_tau1_elecMT", self.MT_tau1_ele)
        self.out.fillBranch("SVFit_tau1_muMT",   self.MT_tau1_mu)
        self.out.fillBranch("SVFit_tau1_hadMT",  self.MT_tau1_tau)
        self.out.fillBranch("SVFit_tau2_muMT",   self.MT_tau2_mu)
        self.out.fillBranch("SVFit_tau2_hadMT",  self.MT_tau2_tau)
        self.out.fillBranch("SVFit_ditauMass", self.tausMass)
        self.out.fillBranch("SVFit_ditauPt", self.ditauPt)
        self.out.fillBranch("SVFit_ditauDR", self.taudr)
        self.out.fillBranch("SVFit_deltaREMu", self.drEMu)

        count = 0                
        for bj in bjet:
            if count == 2: break
            self.out.fillBranch("SVFit_bj" + str(count + 1) + "Pt", bj.pt)
            self.out.fillBranch("SVFit_bj" + str(count + 1) + "Eta", bj.eta)
            self.out.fillBranch("SVFit_bj" + str(count + 1) + "Phi", bj.phi)
            self.out.fillBranch("SVFit_bj" + str(count + 1) + "Mass", bj.mass)
            self.out.fillBranch("SVFit_bj" + str(count + 1) + "DeepCSV", bj.btagDeepB)
            count += 1

        count = 0
        for j in jet:
            if count == 3: break
            if count == 0: jet1.SetPtEtaPhiM(j.pt, j.eta, j.phi, j.mass)
            if count == 1: jet1.SetPtEtaPhiM(j.pt, j.eta, j.phi, j.mass)
            self.out.fillBranch("SVFit_j" + str(count + 1) + "Pt", j.pt)
            self.out.fillBranch("SVFit_j" + str(count + 1) + "Eta", j.eta)
            self.out.fillBranch("SVFit_j" + str(count + 1) + "Phi", j.phi)
            self.out.fillBranch("SVFit_j" + str(count + 1) + "Mass", j.mass)
            self.out.fillBranch("SVFit_j" + str(count + 1) + "DeepCSV", j.btagDeepB)
            count += 1

        self.out.fillBranch("SVFit_2tau2jetPt", self.SVFit_2tau2jetPt)
        self.out.fillBranch("SVFit_dijetPt", addFourVec(jet).Pt())
        self.out.fillBranch("SVFit_dijetMass", addFourVec(jet).M())
        self.out.fillBranch("SVFit_dibjetDR", self.DeltaR(bjet))
        self.out.fillBranch("SVFit_MaxbjetTauDR", self.MaxdrBJetTaus)
        self.out.fillBranch("SVFit_MinbjetTauDR", self.MindrBJetTaus)
        self.out.fillBranch("SVFit_dijetDR", self.DeltaR(jet))
        self.out.fillBranch("SVFit_dijetDEta", self.DeltaEta(jet))
        
        	
        return True
        
        
 # define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed
