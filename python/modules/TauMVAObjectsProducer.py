import os, sys
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import math
import numba
import numpy as np
from array import array
from importlib import import_module

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.tools import deltaPhi, deltaR, closest
from PhysicsTools.NanoAODTools.postprocessing.framework.treeReaderArrayTools import *
from rootpy.tree import Tree, TreeModel, IntCol, FloatArrayCol

from TauAnalysis import ClassicSVfit

@numba.jit(nopython=True)
def recursiveMotherSearch(startIdx, targetIdx, GenPartCut_genPartIdxMother):
    if startIdx < 0:
        return False

    mom = GenPartCut_genPartIdxMother[startIdx]

    if mom < 0:
        return False
    elif startIdx == targetIdx:
        return True
    else:
        return recursiveMotherSearch(mom, targetIdx, GenPartCut_genPartIdxMother)

@numba.jit(nopython=True)
def genParticleAssociation(GenPart_genPartIdxMother, GenPart_pdgId, GenPart_statusFlags):
    GenPart_momPdgId = GenPart_pdgId[GenPart_genPartIdxMother]
    #set any particles with an index of -1 to id 0                                                                                                                      
    GenPart_momPdgId[GenPart_genPartIdxMother < 0] = 0

    genTauDaughters = []
    whichTau = []
    for iGP, pdgId in enumerate(GenPart_pdgId):
        if (GenPart_statusFlags[iGP] & 0x2100) != 0x2100:
            continue
        #print("Mother index {0} and pdgId {1}".format(iGP, pdgId))
        if abs(pdgId) == 15:
            gtd = []
            for iGP2, pdgId2 in enumerate(GenPart_pdgId):
                if abs(pdgId2) >= 11 and abs(pdgId2) <= 14 and (GenPart_statusFlags[iGP2] & 0x2008) == 0x2008:
                    if recursiveMotherSearch(iGP2, iGP, GenPart_genPartIdxMother):
                        if pdgId == 15:    whichTau.append(0)
                        elif pdgId == -15: whichTau.append(1)
                        else:              whichTau.append(-1)
                        gtd.append(iGP2)
	        if (abs(pdgId2) == 211 or abs(pdgId2) == 111) and (GenPart_statusFlags[iGP2] & 0x2008) == 0x2008:
                    if recursiveMotherSearch(iGP2, iGP, GenPart_genPartIdxMother):
                        if pdgId == 15:    whichTau.append(0)
                        elif pdgId == -15: whichTau.append(1)
                        else:              whichTau.append(-1)
                        gtd.append(iGP2)
                #print("Appended index {0} and pdgId {1} and gtd {2}".format(iGP2, pdgId2, gtd))
            genTauDaughters.extend(gtd)

    #print("How many tau daughters? {0}".format(genTauDaughters))
    return genTauDaughters, whichTau

def deltaRMatch(PFCandEta, PFCandPhi, genTauDaughters_eta, genTauDaughters_phi):

    matches = np.zeros(len(PFCandEta), dtype=int)

    if len(genTauDaughters_eta):

        topEtaVals = np.array(np.meshgrid(PFCandEta, genTauDaughters_eta)).T.reshape(-1,2)
        topPhiVals = np.array(np.meshgrid(PFCandPhi, genTauDaughters_phi)).T.reshape(-1,2)
    
        ## Using ufunc for vector operation
        deta = np.power(topEtaVals[:,0] - topEtaVals[:,1], 2)
        dPhi = topPhiVals[:,0] - topPhiVals[:,1]
        #print("DEta {0} and DPhi {1}".format(deta, dPhi))
        dR = np.sqrt((( abs(abs(dPhi)-np.pi)-np.pi )**2+(deta)**2)).reshape([-1,len(genTauDaughters_eta)/len(genTauDaughters_eta), len(genTauDaughters_eta)])

        #print(dR)

        matches[dR.max(axis=2).min(axis=1) < 0.6] = 1
    
    return matches

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
        self.out.branch("nJets30",       "I")
        self.out.branch("HiggsSVFit_PassBaseline", "F")
        self.out.branch("HiggsSVFit_Pt"  ,   "F")
        self.out.branch("HiggsSVFit_Eta" ,   "F")
        self.out.branch("HiggsSVFit_Phi" ,   "F")
        self.out.branch("HiggsSVFit_Mass",   "F")
        self.out.branch("HiggsSVFit_TransverseMass", "F")
        self.out.branch("HiggsSVFit_METRho", "F")
        self.out.branch("HiggsSVFit_METPhi", "F")
        self.out.branch("HiggsSVFit_channel", "F")
        self.out.branch("HiggsSVFit_tau1DM", "I")
        self.out.branch("HiggsSVFit_tau1Mass", "F")
        self.out.branch("HiggsSVFit_tau1pdgId", "I")
        self.out.branch("HiggsSVFit_tau1Pt", "F")
        self.out.branch("HiggsSVFit_tau1Eta", "F")
        self.out.branch("HiggsSVFit_tau1Phi", "F")
        self.out.branch("HiggsSVFit_tau2DM", "I")
        self.out.branch("HiggsSVFit_tau2Mass", "F")
        self.out.branch("HiggsSVFit_tau2pdgId", "I")
        self.out.branch("HiggsSVFit_tau2Pt", "F")
        self.out.branch("HiggsSVFit_tau2Eta", "F")
        self.out.branch("HiggsSVFit_tau2Phi", "F")
        self.out.branch("HiggsSVFit_elecMuonMT",  "F")
        self.out.branch("HiggsSVFit_MaxMT",       "F")
        self.out.branch("HiggsSVFit_tau1_elecMT", "F")
        self.out.branch("HiggsSVFit_tau1_muMT",   "F")
        self.out.branch("HiggsSVFit_tau1_hadMT",  "F")
        self.out.branch("HiggsSVFit_tau2_muMT",   "F")
        self.out.branch("HiggsSVFit_tau2_hadMT",  "F")
        self.out.branch("HiggsSVFit_ditauMass", "F")
        self.out.branch("HiggsSVFit_ditauPt",   "F")
        self.out.branch("HiggsSVFit_ditauDR",    "F")
        self.out.branch("HiggsSVFit_deltaREMu", "F")
        self.out.branch("HiggsSVFit_bj1Pt",   "F")
        self.out.branch("HiggsSVFit_bj1Eta",   "F")
        self.out.branch("HiggsSVFit_bj1Phi",   "F")
        self.out.branch("HiggsSVFit_bj1Mass",   "F")
        self.out.branch("HiggsSVFit_bj1DeepCSV",   "F")
        self.out.branch("HiggsSVFit_bj2Pt",   "F")
        self.out.branch("HiggsSVFit_bj2Eta",   "F")
        self.out.branch("HiggsSVFit_bj2Phi",   "F")
        self.out.branch("HiggsSVFit_bj2Mass",   "F")
        self.out.branch("HiggsSVFit_bj2DeepCSV",   "F")
        self.out.branch("HiggsSVFit_j1Pt",   "F")
        self.out.branch("HiggsSVFit_j1Eta",   "F")
        self.out.branch("HiggsSVFit_j1Phi",   "F")
        self.out.branch("HiggsSVFit_j1Mass",   "F")
        self.out.branch("HiggsSVFit_j1DeepCSV",   "F")
        self.out.branch("HiggsSVFit_j2Pt",   "F")
        self.out.branch("HiggsSVFit_j2Eta",   "F")
        self.out.branch("HiggsSVFit_j2Phi",   "F")
        self.out.branch("HiggsSVFit_j2Mass",   "F")
        self.out.branch("HiggsSVFit_j2DeepCSV",   "F")
        self.out.branch("HiggsSVFit_j3Pt",   "F")
        self.out.branch("HiggsSVFit_j3Eta",   "F")
        self.out.branch("HiggsSVFit_j3Phi",   "F")
        self.out.branch("HiggsSVFit_j3Mass",   "F")
        self.out.branch("HiggsSVFit_j3DeepCSV",   "F")
        self.out.branch("HiggsSVFit_2tau2jetPt", "F")
        self.out.branch("HiggsSVFit_dijetPt",   "F")
        self.out.branch("HiggsSVFit_dijetMass", "F")
        self.out.branch("HiggsSVFit_dijetDR",    "F")
        self.out.branch("HiggsSVFit_dijetDEta",  "F")
        self.out.branch("HiggsSVFit_HT",        "F")


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
        if jet.pt < 30 or math.fabs(jet.eta) > 4.7:
            return False
        return True

    def CalMtWTau1(self, obj, met):
        if obj.tau1Pt != -999: return math.sqrt( 2 * met.pt * obj.tau1Pt * (1 - math.cos(ROOT.TVector2.Phi_mpi_pi(met.phi-obj.tau1Phi))))
        else: return -999

    def CalMtWTau2(self, obj, met):
        if obj.tau2Pt != -999: return math.sqrt( 2 * met.pt * obj.tau2Pt * (1 - math.cos(ROOT.TVector2.Phi_mpi_pi(met.phi-obj.tau2Phi))))
        else: return -999

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

    class TTreeReaderArrayWrapper:
        def __init__(self, ttarray):
            self.ttarray = ttarray

        def __iter__(self):
            for i in xrange(len(self.ttarray)):
                yield self.ttarray[i]
            return

    def nGenParts(self, event):
        GenPart_pdgId = np.fromiter(self.TTreeReaderArrayWrapper(event.GenPart_pdgId), int)
        GenPart_statusFlags = np.fromiter(self.TTreeReaderArrayWrapper(event.GenPart_statusFlags), int)
        GenPart_eta = np.fromiter(self.TTreeReaderArrayWrapper(event.GenPart_eta), float)
        GenPart_phi = np.fromiter(self.TTreeReaderArrayWrapper(event.GenPart_phi), float)

        PFCand_eta = np.fromiter(self.TTreeReaderArrayWrapper(event.PFCands_eta), float)
        PFCand_phi = np.fromiter(self.TTreeReaderArrayWrapper(event.PFCands_phi), float)

        # statusFlag 0x2100 corresponds to "isLastCopy and fromHardProcess"
        # statusFlag 0x2080 corresponds to "IsLastCopy and isHardProcess"
        genPartsFilt = (((abs(GenPart_pdgId) >= 11) & (abs(GenPart_pdgId) <= 14)) | (abs(GenPart_pdgId) == 211) | (abs(GenPart_pdgId) == 111)) & (((GenPart_statusFlags & 0x2008) == 0x2008) | ((GenPart_statusFlags & 0x2004) == 0x2004))

        #calculate deltaR matches
        etas = np.array(np.meshgrid(PFCand_eta, GenPart_eta[genPartsFilt])).T
        deta = np.power(etas[:, :, 0] - etas[:, :, 1], 2)
        phis = np.array(np.meshgrid(PFCand_phi, GenPart_phi[genPartsFilt])).T
        dPhi = phis[:, :, 0] - phis[:, :, 1]
        np.subtract(dPhi, 2*math.pi, out = dPhi, where= (dPhi >=math.pi))
        np.add(dPhi, 2*math.pi,  out =dPhi , where=(dPhi < -1*math.pi))
        np.power(dPhi, 2, out=dPhi)
        dR2 = np.add(deta, dPhi)
        nGenPartMatch = (dR2 < 0.6*0.6).sum(axis=1)
        return nGenPartMatch

    def PFCandGenMatch(self, event, PFCandEta, PFCandPhi):
        
        GenPart_eta              = np.fromiter(self.TTreeReaderArrayWrapper(event.GenPart_eta),              dtype=float)
        GenPart_phi              = np.fromiter(self.TTreeReaderArrayWrapper(event.GenPart_phi),              dtype=float)

        GenPart_genPartIdxMother = np.fromiter(self.TTreeReaderArrayWrapper(event.GenPart_genPartIdxMother), dtype=int)
        GenPart_pdgId            = np.fromiter(self.TTreeReaderArrayWrapper(event.GenPart_pdgId),            dtype=int)
        GenPart_statusFlags      = np.fromiter(self.TTreeReaderArrayWrapper(event.GenPart_statusFlags),      dtype=int)
    
        genTauDaughters_list, whichTaus_list = genParticleAssociation(GenPart_genPartIdxMother, GenPart_pdgId, GenPart_statusFlags)
    
        genTauDaughters, whichTau = np.array(genTauDaughters_list), np.array(whichTaus_list)

    
        if(len(genTauDaughters)):
            genTauDaughters_eta = GenPart_eta[genTauDaughters]
            genTauDaughters_phi = GenPart_phi[genTauDaughters]
        else:
            genTauDaughters_eta = np.array([])
            genTauDaughters_phi = np.array([])

    
        return deltaRMatch(PFCandEta, PFCandPhi, genTauDaughters_eta, genTauDaughters_phi), whichTau

    def SelTaus(self, obj):
        if obj.channel == 0 and (math.fabs(obj.tau1Eta) > 2.1 or obj.tau1Pt < 23 or math.fabs(obj.tau2Eta) > 2.3 or obj.tau2Pt < 30):
            return False
        if obj.channel == 1 and (math.fabs(obj.tau1Eta) > 2.1 or obj.tau1Pt < 26 or math.fabs(obj.tau2Eta) > 2.3 or obj.tau2Pt < 30):
            return False
        if obj.channel == 2 and (math.fabs(obj.tau1Eta) > 2.1 or obj.tau1Pt < 40 or math.fabs(obj.tau2Eta) > 2.1 or obj.tau2Pt < 40):
            return False
        if obj.channel == 3 and (math.fabs(obj.tau1Eta) > 2.4 or obj.tau1Pt < 10 or math.fabs(obj.tau2Eta) > 2.4 or obj.tau2Pt < 10):
            return False
        if obj.channel == 4 and (math.fabs(obj.tau1Eta) > 2.5 or obj.tau1Pt < 13 or math.fabs(obj.tau2Eta) > 2.5 or obj.tau2Pt < 13):
            return False
        if obj.channel == 5 and (math.fabs(obj.tau1Eta) > 2.5 or obj.tau1Pt < 13 or math.fabs(obj.tau2Eta) > 2.4 or obj.tau2Pt < 10):
            return False
        return True

    def analyze(self, event):
        ## Getting objects
        met	  = Object(event, self.metBranchName)
        jets	  = Collection(event, "Jet")
        pfcand    = Collection(event, "PFCands")
        svfit     = Collection(event, "SVFit")
        svfitmet  = Collection(event, "SVFitMET")
        eventNum  = event.event
        
        PFCandPt = np.fromiter(self.TTreeReaderArrayWrapper(event.PFCands_pt), dtype=float)
        PFCandEta = np.fromiter(self.TTreeReaderArrayWrapper(event.PFCands_eta), dtype=float)
        PFCandPhi = np.fromiter(self.TTreeReaderArrayWrapper(event.PFCands_phi), dtype=float)
        
        PFCandGenMatch, whichTau = self.PFCandGenMatch(event, PFCandEta, PFCandPhi)
        
        #print(PFCandGenMatch)
        self.SVFit_Stop0l    = map(lambda x : self.SelTaus(x), svfit)
        self.SVFit_tau1_MtW  = map(lambda x : self.CalMtWTau1(x, met), svfit)
        self.SVFit_tau2_MtW  = map(lambda x : self.CalMtWTau2(x, met), svfit)

        SVIndex = 0
        SVPass = -1
        checkPt = -99
        for sv in xrange(len(svfit)):
            if self.SVFit_Stop0l[sv] and svfit[sv].Pt > checkPt:
                SVPass = sv
                checkPt = svfit[sv].Pt
            elif svfit[sv].Pt > checkPt: 
                SVIndex = sv
                checkPt = svfit[sv].Pt
        SVIndex = SVPass if SVPass != -1 else SVIndex
        #print("Picked Index: {7}, channel: {6}, (pt, eta, phi, mass) = ({0}, {1}, {2}, {3}), Tau1 Pt = {4}, Tau2 Pt = {5}".format(svfit[SVIndex].Pt, svfit[SVIndex].Eta, svfit[SVIndex].Phi, svfit[SVIndex].Mass, svfit[SVIndex].tau1Pt, svfit[SVIndex].tau2Pt, svfit[SVIndex].channel, SVIndex))
                

        
        self.out.fillBranch("HiggsSVFit_PassBaseline",   self.SVFit_Stop0l[SVIndex])
        self.out.fillBranch("HiggsSVFit_Pt"  ,   svfit[SVIndex].Pt)
        self.out.fillBranch("HiggsSVFit_Eta" ,   svfit[SVIndex].Eta)
        self.out.fillBranch("HiggsSVFit_Phi" ,   svfit[SVIndex].Phi)
        self.out.fillBranch("HiggsSVFit_Mass",   svfit[SVIndex].Mass)
        self.out.fillBranch("HiggsSVFit_TransverseMass",   svfit[SVIndex].TransverseMass)
        self.out.fillBranch("HiggsSVFit_METRho", svfit[SVIndex].METRho)
        self.out.fillBranch("HiggsSVFit_METPhi", svfit[SVIndex].METPhi)
        self.out.fillBranch("HiggsSVFit_channel", svfit[SVIndex].channel)
        self.out.fillBranch("HiggsSVFit_tau1DM", int(svfit[SVIndex].tau1DM))
        self.out.fillBranch("HiggsSVFit_tau1pdgId", int(svfit[SVIndex].tau1pdgId))
        self.out.fillBranch("HiggsSVFit_tau1Mass", svfit[SVIndex].tau1Mass)
        self.out.fillBranch("HiggsSVFit_tau1Pt", svfit[SVIndex].tau1Pt)
        self.out.fillBranch("HiggsSVFit_tau1Eta", svfit[SVIndex].tau1Eta)
        self.out.fillBranch("HiggsSVFit_tau1Phi", svfit[SVIndex].tau1Phi)
        self.out.fillBranch("HiggsSVFit_tau2DM", int(svfit[SVIndex].tau2DM))
        self.out.fillBranch("HiggsSVFit_tau2pdgId", int(svfit[SVIndex].tau2pdgId))
        self.out.fillBranch("HiggsSVFit_tau2Mass", svfit[SVIndex].tau2Mass)
        self.out.fillBranch("HiggsSVFit_tau2Pt", svfit[SVIndex].tau2Pt)
        self.out.fillBranch("HiggsSVFit_tau2Eta", svfit[SVIndex].tau2Eta)
        self.out.fillBranch("HiggsSVFit_tau2Phi", svfit[SVIndex].tau2Phi)

        drEMu = deltaR(svfit[SVIndex].tau1Eta, svfit[SVIndex].tau1Phi, svfit[SVIndex].tau2Eta, svfit[SVIndex].tau2Phi) if svfit[SVIndex].channel == 5 else -999
        taus = ROOT.TLorentzVector()
        tau1 = ROOT.TLorentzVector()
        tau2 = ROOT.TLorentzVector()
        MET  = ROOT.TLorentzVector()
        tau1.SetPtEtaPhiM(svfit[SVIndex].tau1Pt, svfit[SVIndex].tau1Eta, svfit[SVIndex].tau1Phi, svfit[SVIndex].tau1Mass)
        tau2.SetPtEtaPhiM(svfit[SVIndex].tau2Pt, svfit[SVIndex].tau2Eta, svfit[SVIndex].tau2Phi, svfit[SVIndex].tau2Mass)
        taus = (tau1 + tau2)
        taudr = deltaR(svfit[SVIndex].tau1Eta, svfit[SVIndex].tau1Phi, svfit[SVIndex].tau2Eta, svfit[SVIndex].tau2Phi)
        MT_elecMu   = self.CalMtW(taus, met) if svfit[SVIndex].channel == 5 else -999
        MT_tau1_ele = self.SVFit_tau1_MtW[SVIndex] if svfit[SVIndex].channel == 5 or svfit[SVIndex].channel == 1 else -999
        MT_tau1_mu  = self.SVFit_tau1_MtW[SVIndex] if svfit[SVIndex].channel == 0 or svfit[SVIndex].channel == 3 else -999
        MT_tau1_tau = self.SVFit_tau1_MtW[SVIndex] if svfit[SVIndex].channel == 2 else -999
        MT_tau2_mu  = self.SVFit_tau2_MtW[SVIndex] if (svfit[SVIndex].channel == 5 or svfit[SVIndex].channel == 3) else -999
        MT_tau2_tau = self.SVFit_tau2_MtW[SVIndex] if svfit[SVIndex].channel <= 2 else -999

        self.out.fillBranch("HiggsSVFit_elecMuonMT", MT_elecMu)
        self.out.fillBranch("HiggsSVFit_MaxMT", max(MT_tau1_ele, MT_tau2_mu))
        self.out.fillBranch("HiggsSVFit_tau1_elecMT", MT_tau1_ele)
        self.out.fillBranch("HiggsSVFit_tau1_muMT",   MT_tau1_mu)
        self.out.fillBranch("HiggsSVFit_tau1_hadMT",  MT_tau1_tau)
        self.out.fillBranch("HiggsSVFit_tau2_muMT",   MT_tau2_mu)
        self.out.fillBranch("HiggsSVFit_tau2_hadMT",  MT_tau2_tau)
        self.out.fillBranch("HiggsSVFit_ditauMass", taus.M())
        MET.SetPtEtaPhiM(met.pt, 0., met.phi, 0.)
        taus = (tau1 + tau2 + MET)
        self.out.fillBranch("HiggsSVFit_ditauPt", taus.Pt())
        self.out.fillBranch("HiggsSVFit_ditauDR", taudr)
        self.out.fillBranch("HiggsSVFit_deltaREMu", drEMu)

        #Jet Variables
        self.Jets_Stop0l      = map(self.SelJets, jets)
        self.out.fillBranch("nJets30",          sum(self.Jets_Stop0l))
        HT = self.CalHT(jets)
        self.out.fillBranch("HiggsSVFit_HT",       HT)
        jetvec = ROOT.TLorentzVector()
        jet1 = ROOT.TLorentzVector()
        jet2 = ROOT.TLorentzVector()
        jet = []
        bjet = []
        for j in xrange(len(jets)):
            if self.Jets_Stop0l[j]: jet.append(jets[j])
            if jets[j].btagStop0l: bjet.append(jets[j])

        count = 0                
        for bj in bjet:
            if count == 2: break
            self.out.fillBranch("HiggsSVFit_bj" + str(count + 1) + "Pt", bj.pt)
            self.out.fillBranch("HiggsSVFit_bj" + str(count + 1) + "Eta", bj.eta)
            self.out.fillBranch("HiggsSVFit_bj" + str(count + 1) + "Phi", bj.phi)
            self.out.fillBranch("HiggsSVFit_bj" + str(count + 1) + "Mass", bj.mass)
            self.out.fillBranch("HiggsSVFit_bj" + str(count + 1) + "DeepCSV", bj.btagDeepB)
            count += 1

        count = 0
        for j in jet:
            if count == 3: break
            if count == 0: jet1.SetPtEtaPhiM(j.pt, j.eta, j.phi, j.mass)
            if count == 1: jet1.SetPtEtaPhiM(j.pt, j.eta, j.phi, j.mass)
            self.out.fillBranch("HiggsSVFit_j" + str(count + 1) + "Pt", j.pt)
            self.out.fillBranch("HiggsSVFit_j" + str(count + 1) + "Eta", j.eta)
            self.out.fillBranch("HiggsSVFit_j" + str(count + 1) + "Phi", j.phi)
            self.out.fillBranch("HiggsSVFit_j" + str(count + 1) + "Mass", j.mass)
            self.out.fillBranch("HiggsSVFit_j" + str(count + 1) + "DeepCSV", bj.btagDeepB)
            count += 1

        jetvec = (jet1 + jet2 + tau1 + tau2 + MET)
        self.out.fillBranch("HiggsSVFit_2tau2jetPt", jetvec.Pt())
        self.out.fillBranch("HiggsSVFit_dijetPt", self.addFourVec(jet).Pt())
        self.out.fillBranch("HiggsSVFit_dijetMass", self.addFourVec(jet).M())
        self.out.fillBranch("HiggsSVFit_dijetDR", self.DeltaR(jet))
        self.out.fillBranch("HiggsSVFit_dijetDEta", self.DeltaEta(jet))
        
        	
        return True
        
        
 # define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed
