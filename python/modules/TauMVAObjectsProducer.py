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

#2016 MC: https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation80XReReco#Data_MC_Scale_Factors_period_dep
#2017 MC: https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation94X

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
        self.out.branch("Tau_channel",   "I")
        self.out.branch("nTauMatch",     "I")
        self.out.branch("Tau_hadDecayFlag1", "I")
        self.out.branch("Tau_hadDecayFlag2", "I")
        self.out.branch("Jet_dijetMass", "F")
        self.out.branch("Tau_dijetMass", "F")
        self.out.branch("Jet_deltaR",    "F")
        self.out.branch("Tau_deltaR",    "F")
        self.out.branch("HiggsCand_pt",  "F") 
        self.out.branch("HiggsCand_eta", "F") 
        self.out.branch("HiggsCand_phi", "F") 
        self.out.branch("HiggsCand_mass","F") 
        self.out.branch("Tau_HT",        "F")
        self.out.branch("Jet_matchPt_1"  , "F")
        self.out.branch("Jet_matchEta_1" , "F")
        self.out.branch("Jet_matchPhi_1" , "F")
        self.out.branch("Jet_matchMass_1", "F")
        self.out.branch("Jet_matchPt_2"  , "F")
        self.out.branch("Jet_matchEta_2" , "F")
        self.out.branch("Jet_matchPhi_2" , "F")
        self.out.branch("Jet_matchMass_2", "F")
        self.out.branch("Jet_matchPt_3"  , "F")
        self.out.branch("Jet_matchEta_3" , "F")
        self.out.branch("Jet_matchPhi_3" , "F")
        self.out.branch("Jet_matchMass_3", "F")
        self.out.branch("Tau_nLep",  "I")
        self.out.branch("Tau_LeptonPt_1", "F")
        self.out.branch("Tau_LeptonEta_1", "F")
        self.out.branch("Tau_LeptonPhi_1", "F")
        self.out.branch("Tau_LeptonCharge_1", "I")
        self.out.branch("Tau_LeptonPdgId_1", "I")
        self.out.branch("Tau_LeptonPt_2", "F")
        self.out.branch("Tau_LeptonEta_2", "F")
        self.out.branch("Tau_LeptonPhi_2", "F")
        self.out.branch("Tau_LeptonCharge_2", "I")
        self.out.branch("Tau_LeptonPdgId_2", "I")


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

    def CalHT(self, jets):
        HT = sum([j.pt for i, j in enumerate(jets) if self.Jet_Stop0l[i]])
        return HT

    def DeltaR(self, obj):
        if len(obj) > 1: dr = deltaR(obj[0].eta, obj[0].phi, obj[1].eta, obj[1].phi)
        elif len(obj) == 1: dr = deltaR(obj[0].eta, obj[0].phi, 0, 0)
        else: dr = -1
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
    
        genTauDaughters_list, whichTau_list = genParticleAssociation(GenPart_genPartIdxMother, GenPart_pdgId, GenPart_statusFlags)
    
        genTauDaughters, whichTau = np.array(genTauDaughters_list), np.array(whichTau_list)

    
        if(len(genTauDaughters)):
            genTauDaughters_eta = GenPart_eta[genTauDaughters]
            genTauDaughters_phi = GenPart_phi[genTauDaughters]
        else:
            genTauDaughters_eta = np.array([])
            genTauDaughters_phi = np.array([])

    
        return deltaRMatch(PFCandEta, PFCandPhi, genTauDaughters_eta, genTauDaughters_phi), whichTau

    def analyze(self, event):
        ## Getting objects
	met	  = Object(event, self.metBranchName)
	jets	  = Collection(event, "Jet")
	pfcand    = Collection(event, "PFCands")
	eventNum  = event.event

        PFCandPt = np.fromiter(self.TTreeReaderArrayWrapper(event.PFCands_pt), dtype=float)
        PFCandEta = np.fromiter(self.TTreeReaderArrayWrapper(event.PFCands_eta), dtype=float)
        PFCandPhi = np.fromiter(self.TTreeReaderArrayWrapper(event.PFCands_phi), dtype=float)

        PFCandGenMatch, whichTau = self.PFCandGenMatch(event, PFCandEta, PFCandPhi)

        #print(PFCandGenMatch)

        vec = [ROOT.TLorentzVector(), ROOT.TLorentzVector()]
        decayType = [-1, -1]
        hadDecayType = [0, 0]
        hvec = ROOT.TLorentzVector()
        #print(PFCandGenMatch)
        #print(whichTau)
        lep = []
        for pf, pI, mTau in zip(pfcand, PFCandGenMatch, whichTau):
            trk = ROOT.TLorentzVector()
            if abs(pf.pdgId) == 11 or abs(pf.pdgId) == 13:
                trk.SetPtEtaPhiM(pf.pt, pf.eta, pf.phi, pf.mass)
                lep.append(pf)
                if mTau >= 0: decayType[mTau] = 0
	    if (abs(pf.pdgId) == 211 or abs(pf.pdgId) == 111):
                trk.SetPtEtaPhiM(pf.pt, pf.eta, pf.phi, pf.mass)
                if mTau >= 0: 
                    decayType[mTau] = 1
                    if abs(pf.pdgId) == 211: hadDecayType[mTau] += 1

            vec[mTau] += (trk)

        hvec = (vec[0] + vec[1])
        #print("Tau 1 ({0}, {1}, {2}, {3}) and Tau 2 ({4}, {5}, {6}, {7})".format(vec[0].Pt(), vec[0].Eta(), vec[0].Phi(), vec[0].M(), vec[1].Pt(), vec[1].Eta(), vec[1].Phi(), vec[1].M()))

        taudr = deltaR(vec[0].Eta(), vec[0].Phi(), vec[1].Eta(), vec[1].Phi()) if sum(PFCandGenMatch) > 0 else -9
        taumass = hvec if sum(PFCandGenMatch) > 0 else ROOT.TLorentzVector(-9, 0, 0, -9)
        self.out.fillBranch("Tau_channel", sum(decayType))
        self.out.fillBranch("nTauMatch",     len(hadDecayType))
        self.out.fillBranch("Tau_hadDecayFlag1", hadDecayType[0])
        self.out.fillBranch("Tau_hadDecayFlag2", hadDecayType[1])
        self.out.fillBranch("Tau_dijetMass", taumass.M())
        self.out.fillBranch("Tau_deltaR", taudr)

        #Lepton Info
        self.out.fillBranch("Tau_nLep", len(lep))
        for l in xrange(len(lep)):
            if l == 2: break
            self.out.fillBranch("Tau_LeptonPt_" + str(l + 1), lep[l].pt)
            self.out.fillBranch("Tau_LeptonEta_" + str(l + 1), lep[l].eta)
            self.out.fillBranch("Tau_LeptonPhi_" + str(l + 1), lep[l].phi)
            self.out.fillBranch("Tau_LeptonCharge_" + str(l + 1), lep[l].charge)
            self.out.fillBranch("Tau_LeptonPdgId_" + str(l + 1), lep[l].pdgId)

        #Add met to di-tau 4-vector
        metvec = ROOT.TLorentzVector()
        metvec.SetPtEtaPhiM(met.pt, 0, met.phi, 0)
        hvec += (metvec)

        self.out.fillBranch("HiggsCand_pt",   hvec.Pt())
        self.out.fillBranch("HiggsCand_eta",  hvec.Eta())
        self.out.fillBranch("HiggsCand_phi",  hvec.Phi())
        self.out.fillBranch("HiggsCand_mass", hvec.M())

        #self.out.fillBranch("nGenHadTaus", 	nGenHadTaus)
	#self.out.fillBranch("nGenTaus", 	nGenTaus)
	#self.out.fillBranch("nGenChHads", 	nGenChHads)
	#self.out.fillBranch("nGenLeptons", 	nGenLeptons)

        #Jet Variables
        self.Jet_Stop0l      = map(self.SelJets, jets)
        self.out.fillBranch("nJets30",          sum(self.Jet_Stop0l))
        HT = self.CalHT(jets)
        self.out.fillBranch("Tau_HT",       HT)
        jetvec = ROOT.TLorentzVector()
        jet = []
        jetpt = []
        jeteta = []
        jetphi = []
        jetmass = []
        for j in xrange(len(jets)):
            if j == 3: break
            if self.Jet_Stop0l[j]:
                jet.append(jets[j])
                self.out.fillBranch("Jet_matchPt_"  + str(j + 1), jets[j].pt)
                self.out.fillBranch("Jet_matchEta_" + str(j + 1), jets[j].eta)
                self.out.fillBranch("Jet_matchPhi_" + str(j + 1), jets[j].phi)
                self.out.fillBranch("Jet_matchMass_"+ str(j + 1), jets[j].mass)

        self.out.fillBranch("Jet_dijetMass", self.addFourVec(jet).M())
        self.out.fillBranch("Jet_deltaR", self.DeltaR(jet))

		
	return True


 # define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed
