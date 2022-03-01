import os, sys
from os import path
import io
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import math
import numba
import numpy as np
import gzip
from array import array
from importlib import import_module

from PhysicsTools.NanoHiggsTools.modules.HelperFunctions import *
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.tools import deltaPhi, deltaR, closest
from PhysicsTools.NanoAODTools.postprocessing.framework.treeReaderArrayTools import *
from rootpy.tree import Tree, TreeModel, IntCol, FloatArrayCol

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
    for iGP, pdgId in enumerate(GenPart_pdgId):
        if (GenPart_statusFlags[iGP] & 0x2100) != 0x2100:
            continue
        if abs(pdgId) == 15:
            gtd = []
            for iGP2, pdgId2 in enumerate(GenPart_pdgId):
                if (abs(pdgId2) == 11 or abs(pdgId2) == 13 or abs(pdgId2) == 211) and (GenPart_statusFlags[iGP2] & 0x2008) == 0x2008:
                    if recursiveMotherSearch(iGP2, iGP, GenPart_genPartIdxMother):
                        gtd.append(iGP2)
            genTauDaughters.extend(gtd)

    return genTauDaughters

@numba.jit(nopython=True)
def isotrkAssociation(PFCands_pdgId):
    genTauDaughters = []
    genTaupdgId = []
    gtd = []
    gtpdgid = []
    for iGP, pdgId in enumerate(PFCands_pdgId):
        if (abs(pdgId) == 11 or abs(pdgId) == 13 or abs(pdgId) == 211):
            gtd.append(iGP)
            gtpdgid.append(pdgId)
    genTauDaughters.extend(gtd)
    genTaupdgId.extend(gtpdgid)

    return genTauDaughters, genTaupdgId

def deltaRMatch(CandPdgId, CandEta, CandPhi, genTauDaughters_pdgId, genTauDaughters_eta, genTauDaughters_phi):

    matches = np.zeros(len(CandEta), dtype=int)
    candmatches = np.zeros(len(genTauDaughters_eta), dtype=int)
    matchesPdgId = np.zeros(len(CandEta), dtype=int)

    if len(genTauDaughters_eta):

        EtaVals = np.array(np.meshgrid(CandEta, genTauDaughters_eta)).T.reshape(-1,2)
        PhiVals = np.array(np.meshgrid(CandPhi, genTauDaughters_phi)).T.reshape(-1,2)
        PdgIdVals = np.array(np.meshgrid(CandPdgId, genTauDaughters_pdgId)).T.reshape(-1,2)

        #print("candEta: {0}, genTauEta: {1}".format(CandEta, genTauDaughters_eta))

        ## Using ufunc for vector operation
        deta = np.power(EtaVals[:,0] - EtaVals[:,1], 2)
        dPhi = PhiVals[:,0] - PhiVals[:,1]
        dR = np.sqrt((( abs(abs(dPhi)-np.pi)-np.pi )**2+(deta)**2)).reshape([-1, len(genTauDaughters_eta), 1])
        dRMatch = dR.T.reshape([-1, len(CandEta), 1])
        ## Check PDG agreement
        matchPdgID = (PdgIdVals[:,0] == abs(PdgIdVals[:,1])).reshape([-1, len(genTauDaughters_eta), 1])

        #print("dR: {0}, dRMatches: {1}".format(dR, dRMatch))

        matches[dR.max(axis=2).min(axis=1) < 0.8] = 1
        candmatches[dRMatch.max(axis=2).min(axis=1) < 0.8] = 1
        matchesPdgId[matchPdgID.max(axis=2).min(axis=1)] = 1
        #print("matches: {0}, candmatches: {1}".format(matches, candmatches))

    return matches, candmatches, matchesPdgId

TauDecay = {
    -999. : 'kUnknown',
    -111. : 'kUnknown',
    1 : 'kTauToHadDecay',
    2 : 'kTauToElecDecay',
    3 : 'kTauToMuDecay',
}

tauDecayPdgId = {
    'kUnknown' : -1,
    'kTauToHadDecay' : 211,
    'kTauToElecDecay' : 11,
    'kTauToMuDecay' : 13,
}

class HiggsLundVarsProducer(Module):
    def __init__(self, processName, match):
        self.metBranchName = "MET"
        self.processName = processName
        self.match = match
        self.debug = False

    def beginJob(self):
        pass
    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("SVFit_LundHiggsM",         "F", lenVar="2*nSVFit")
        self.out.branch("SVFit_LundHiggsKT",        "F", lenVar="2*nSVFit")
        self.out.branch("SVFit_LundHiggsZ",         "F", lenVar="2*nSVFit")
        self.out.branch("SVFit_LundHiggsDelta",     "F", lenVar="2*nSVFit")
        self.out.branch("SVFit_LundHiggsKappa",     "F", lenVar="2*nSVFit")
        self.out.branch("SVFit_LundHiggsPsi",       "F", lenVar="2*nSVFit")
        self.out.branch("SVFit_LundTau1M",          "F", lenVar="nSVFit")
        self.out.branch("SVFit_LundTau1KT",         "F", lenVar="nSVFit")
        self.out.branch("SVFit_LundTau1Z",          "F", lenVar="nSVFit")
        self.out.branch("SVFit_LundTau1Delta",      "F", lenVar="nSVFit")
        self.out.branch("SVFit_LundTau1Kappa",      "F", lenVar="nSVFit")
        self.out.branch("SVFit_LundTau1Psi",        "F", lenVar="nSVFit")
        self.out.branch("SVFit_LundTau2M",          "F", lenVar="nSVFit")
        self.out.branch("SVFit_LundTau2KT",         "F", lenVar="nSVFit")
        self.out.branch("SVFit_LundTau2Z",          "F", lenVar="nSVFit")
        self.out.branch("SVFit_LundTau2Delta",      "F", lenVar="nSVFit")
        self.out.branch("SVFit_LundTau2Kappa",      "F", lenVar="nSVFit")
        self.out.branch("SVFit_LundTau2Psi",        "F", lenVar="nSVFit")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass
    def makeTLorentzVector(self, svfit, isTau=""):
        v = ROOT.TLorentzVector()
        if isTau=="tau1":   v.SetPtEtaPhiM(svfit.tau1Pt, svfit.tau1Eta, svfit.tau1Phi, svfit.tau1Mass)
        elif isTau=="tau2": v.SetPtEtaPhiM(svfit.tau2Pt, svfit.tau2Eta, svfit.tau2Phi, svfit.tau2Mass)
        elif isTau=="gen":  v.SetPtEtaPhiM(svfit.pt, svfit.eta, svfit.phi, svfit.mass)
        else:               v.SetPtEtaPhiM(svfit.Pt, svfit.Eta, svfit.Phi, svfit.Mass)
        return v

    class TTreeReaderArrayWrapper:
        def __init__(self, ttarray):
            self.ttarray = ttarray

        def __iter__(self):
            for i in xrange(len(self.ttarray)):
                yield self.ttarray[i]
            return

    def HiggsGenMatch(self, event, CandEta, CandPhi, CandPdgId, match = 'GenPart'):
        genTauDaughters_list = []
        genTaupdgId_list = []
    
        eta              = np.fromiter(self.TTreeReaderArrayWrapper(event.GenPart_eta),              dtype=float)
        phi              = np.fromiter(self.TTreeReaderArrayWrapper(event.GenPart_phi),              dtype=float)
        genPartIdxMother = np.fromiter(self.TTreeReaderArrayWrapper(event.GenPart_genPartIdxMother), dtype=int)
        pdgId            = np.fromiter(self.TTreeReaderArrayWrapper(event.GenPart_pdgId),            dtype=int)
        statusFlags      = np.fromiter(self.TTreeReaderArrayWrapper(event.GenPart_statusFlags),      dtype=int)
    
        genTauDaughters_list = genParticleAssociation(genPartIdxMother, pdgId, statusFlags)
    
        genTauDaughters = np.array(genTauDaughters_list)
        genTaupdgId = np.array(genTaupdgId_list)
    
        if(len(genTauDaughters)):
            genTauDaughters_pdgid = pdgId[genTauDaughters]
            genTauDaughters_eta = eta[genTauDaughters]
            genTauDaughters_phi = phi[genTauDaughters]
        else:
            genTauDaughters_pdgid = np.array([])
            genTauDaughters_eta = np.array([])
            genTauDaughters_phi = np.array([])
    
        matches, candmatches, matchPdgID = deltaRMatch(CandPdgId, CandEta, CandPhi, genTauDaughters_pdgid, genTauDaughters_eta, genTauDaughters_phi)
    
        return matches, candmatches, matchPdgID, genTauDaughters

    def recursiveFindHiggs(self, startIdx, targetPdgId, genpart):
        if startIdx < 0:
            return -1
    
        mom = genpart[startIdx].genPartIdxMother
    
        if mom < 0:
            return startIdx
        elif genpart[startIdx].pdgId == targetPdgId:
            return mom
        else:
            return self.recursiveFindHiggs(mom, targetPdgId, genpart)

    def findwhichTau(self, tau, genpart):
        idx = -1
        PT = -999
        for i, t in enumerate(tau):
            pdgid = genpart[t].pdgId
            pt = genpart[t].pt

            if pt >= PT:
                PT = pt
                idx = i

            #if abs(pdgid) == 25:
            #    return i

        return i

    def whichPT(self, particle, which = "higgs"):
        if which == "tau1":
            return particle.tau1Pt, particle.tau1Eta, particle.tau1Phi
        elif which == "tau2":
            return particle.tau2Pt, particle.tau2Eta, particle.tau2Phi
        else:
            return particle.Pt, particle.Eta, particle.Phi

    def lnKT(self, gen, particle, delta, which = "higgs"):
        pt, eta, phi = self.whichPT(particle, which)

        return np.float32(math.log(pt * delta))

    def lundDelta(self, gen, particle, which = "higgs"):
        pt, eta, phi = self.whichPT(particle, which)

        return np.float32(max(1e-6, deltaR(gen.eta, gen.phi, eta, pt)))

    def Z(self, gen, particle, which = "higgs"):
        pt, eta, phi = self.whichPT(particle, which)
        return np.float32(pt / (gen.pt + pt))

    def rap(self, p):
        MaxRap = 1e5
        if (p.Pt() == 0.0):
          _phi = 0.0
        else:
          _phi = math.atan2(p.Py(),p.Px())
        
        if (_phi < 0.0): _phi += 2*np.pi
        if (_phi >= 2*np.pi): _phi -= 2*np.pi # can happen if phi=-|eps<1e-15|?
        if (p.E() == abs(p.Pz()) and p.Pt() == 0):
          MaxRapHere = MaxRap + abs(p.Pz());
          _rap = MaxRapHere if p.Pz() >= 0.0 else -MaxRapHere
        else:
          effective_m2 = max(0.0,p.M()*p.M())
          E_plus_pz    = p.E() + abs(p.Pz())
          
          _rap = 0.5*math.log((p.Pt() + effective_m2)/(E_plus_pz*E_plus_pz));
          if (p.Pz() > 0):_rap = - _rap

        return _rap


    def Psi(self, gen, particle, which = "higgs"):
        p = self.makeTLorentzVector(particle, which)
        g = self.makeTLorentzVector(gen, "gen")
        rap1 = self.rap(g)
        rap2 = self.rap(p)
        pt, eta, phi = self.whichPT(particle, which)

        try:
            psi = np.float32(math.atan((rap1 - rap2) / (gen.phi - phi)))
        except ZeroDivisionError:
            psi = 0

        return psi

    def lnM(self, gen, particle, which = "higgs"):
        g = self.makeTLorentzVector(gen, "gen")
        p = self.makeTLorentzVector(particle, which)

        m = (g + p).M()

        return np.float32(0.5 * math.log(abs(m*m)))

    def analyze(self, event):
        ## Getting objects
        svfit      = Collection(event, "SVFit")
        svfitmet   = Collection(event, "SVFitMET")
        genpart    = Collection(event, "GenPart")
        genvistau  = Collection(event, "GenVisTau")
        isotrk     = Collection(event, "IsoTrack")
        pfcand     = Collection(event, "PFCands")
        fjtopfcand = Collection(event, "FatJetToPFCands")

        tau1Pt = np.fromiter(self.TTreeReaderArrayWrapper(event.SVFit_tau1Pt), dtype=float)
        tau1Eta = np.fromiter(self.TTreeReaderArrayWrapper(event.SVFit_tau1Eta), dtype=float)
        tau1Phi = np.fromiter(self.TTreeReaderArrayWrapper(event.SVFit_tau1Phi), dtype=float)
        tau1pdgId = np.fromiter(self.TTreeReaderArrayWrapper(event.SVFit_tau1pdgId), dtype=float)
        tau1DM = np.fromiter(self.TTreeReaderArrayWrapper(event.SVFit_tau1DM), dtype=float)
        tau1pdgIdTranslate = [tauDecayPdgId[TauDecay[t]] for t in tau1pdgId]
        tau1GenMatch, tau1GenPartMatch, tau1GenMatchPdgId, whichTau1 = self.HiggsGenMatch(event, tau1Eta, tau1Phi, tau1pdgIdTranslate, self.match)

        tau2Pt = np.fromiter(self.TTreeReaderArrayWrapper(event.SVFit_tau2Pt), dtype=float)
        tau2Eta = np.fromiter(self.TTreeReaderArrayWrapper(event.SVFit_tau2Eta), dtype=float)
        tau2Phi = np.fromiter(self.TTreeReaderArrayWrapper(event.SVFit_tau2Phi), dtype=float)
        tau2pdgId = np.fromiter(self.TTreeReaderArrayWrapper(event.SVFit_tau2pdgId), dtype=float)
        tau2DM = np.fromiter(self.TTreeReaderArrayWrapper(event.SVFit_tau2DM), dtype=float)
        tau2pdgIdTranslate = [tauDecayPdgId[TauDecay[t]] for t in tau2pdgId]
        tau2GenMatch, tau2GenPartMatch, tau2GenMatchPdgId, whichTau2 = self.HiggsGenMatch(event, tau2Eta, tau2Phi, tau2pdgIdTranslate, self.match)

        if self.debug:
            print("tau1 pt: {0} eta: {1} pdgId: {2} DM: {3}".format(tau1Pt, tau1Eta, tau1pdgId, tau1DM))
            print("tau2 pt: {0} eta: {1} pdgId: {2} DM: {3}".format(tau2Pt, tau2Eta, tau2pdgId, tau2DM))
            print("tau1: {0} {1} {2} {3}".format(tau1GenMatch, tau1pdgIdTranslate, whichTau1, tau1GenPartMatch))
            print("tau2: {0} {1} {2} {3}".format(tau2GenMatch, tau2pdgIdTranslate, whichTau2, tau2GenPartMatch))

        tauIdx1 = []
        tauIdx2 = []
        matchtau1 = []
        matchtau2 = []
        higgs1Idx = []
        higgs2Idx = []
        for idx, (gD1, gD2, t1, t2) in enumerate(zip(whichTau1, whichTau2, tau1GenPartMatch, tau2GenPartMatch)):
            idx1 = self.recursiveFindHiggs(gD1, 15, genpart) if t1 else -1
            idx2 = self.recursiveFindHiggs(gD2, 15, genpart) if t2 else -1
            if idx1 >= 0: tauIdx1.append(idx1)
            if idx2 >= 0: tauIdx2.append(idx2)
            if idx1 >= 0:
                higgs1Idx.append(genpart[idx1].genPartIdxMother)
                if self.debug: print("Matched Higgs Index: {0} --> tau1: {1} pdgId: {2}".format(genpart[idx1].genPartIdxMother, idx1, genpart[idx1].pdgId))
                break
            elif idx2 >= 0:
                higgs2Idx.append(genpart[idx2].genPartIdxMother)
                if self.debug: print("Matched Higgs Index: {0} --> tau2: {1} pdgId: {2}".format(genpart[idx2].genPartIdxMother, idx2, genpart[idx2].pdgId))
                break
#---
        if len(tauIdx1) == 0 and len(whichTau1) > 0: 
            tauIdx1.append(genpart[whichTau1[0]].genPartIdxMother)
            higgs1Idx.append(genpart[genpart[whichTau1[0]].genPartIdxMother].genPartIdxMother)
        if len(tauIdx2) == 0 and len(whichTau2) > 0: 
            tauIdx2.append(genpart[whichTau2[0]].genPartIdxMother)
            higgs2Idx.append(genpart[genpart[whichTau2[0]].genPartIdxMother].genPartIdxMother)

        if self.debug: print("Higgs Vec: {0} = {1}, tau1 Vec: {2}, tau2 Vec: {3}".format(higgs1Idx, higgs2Idx, tauIdx1, tauIdx2))

 
        tau1LnM = tau1LnKT = tau1LnZ = tau1LnDelta = tau1LnKappa = tau1LnPsi = []
        tau2LnM = tau2LnKT = tau2LnZ = tau2LnDelta = tau2LnKappa = tau2LnPsi = []
        higgsLnM = higgsLnKT = higgsLnZ = higgsLnDelta = higgsLnKappa = higgsLnPsi = []
        for i, (gD1, gD2, t) in enumerate(zip(tau1GenMatch, tau2GenMatch, svfit)):
            higgs = -1
            if len(higgs1Idx) != 0: higgs = higgs1Idx[0]
            elif len(higgs2Idx) != 0: higgs = higgs2Idx[0]
            if higgs == -1 and len(higgs2Idx) != 0: higgs = higgs2Idx[0]
            if higgs == -1: higgs = 0
            ### Fill for tau1
            if gD1:
                delta = self.lundDelta(genpart[tauIdx1[0]], t, "tau1")
                z = self.Z(genpart[tauIdx1[0]], t, "tau1")
                tau1LnM.append(self.lnM(genpart[tauIdx1[0]], t, "tau1"))
                tau1LnKT.append(self.lnKT(genpart[tauIdx1[0]], t, delta, "tau1"))
                tau1LnZ.append(np.float32(math.log(z)))
                tau1LnDelta.append(np.float32(math.log(delta)))
                tau1LnKappa.append(np.float32(math.log(z * delta)))
                tau1LnPsi.append(self.Psi(genpart[tauIdx1[0]], t, "tau1"))
                ##Fill Higgs for tau 1
                delta = self.lundDelta(genpart[higgs], t, "tau1")
                z = self.Z(genpart[higgs], t, "tau1")
                higgsLnM.append(self.lnM(genpart[higgs], t, "tau1"))
                higgsLnKT.append(self.lnKT(genpart[higgs], t, delta, "tau1"))
                higgsLnZ.append(np.float32(math.log(z)))
                higgsLnDelta.append(np.float32(math.log(delta)))
                higgsLnKappa.append(np.float32(math.log(z * delta)))
                higgsLnPsi.append(self.Psi(genpart[higgs], t, "tau1"))
            else:
                tau1LnM.append(-99999.)
                tau1LnKT.append(-99999.)
                tau1LnZ.append(-99999.)
                tau1LnDelta.append(-99999.)
                tau1LnKappa.append(-99999.)
                tau1LnPsi.append(-99999.)
                higgsLnM.append(-99999.)
                higgsLnKT.append(-99999.)
                higgsLnZ.append(-99999.)
                higgsLnDelta.append(-99999.)
                higgsLnKappa.append(-99999.)
                higgsLnPsi.append(-99999.)
            ### Fill for tau
            if gD2:
                delta = self.lundDelta(genpart[tauIdx1[0]], t, "tau2")
                z = self.Z(genpart[tauIdx1[0]], t, "tau2")
                tau2LnM.append(self.lnM(genpart[tauIdx1[0]], t, "tau2"))
                tau2LnKT.append(self.lnKT(genpart[tauIdx1[0]], t, delta, "tau2"))
                tau2LnZ.append(np.float32(math.log(z)))
                tau2LnDelta.append(np.float32(math.log(delta)))
                tau2LnKappa.append(np.float32(math.log(z * delta)))
                tau2LnPsi.append(self.Psi(genpart[tauIdx1[0]], t, "tau2"))
                ##Fill Higgs for tau 1
                delta = self.lundDelta(genpart[higgs], t, "tau2")
                z = self.Z(genpart[higgs], t, "tau2")
                higgsLnM.append(self.lnM(genpart[higgs], t, "tau2"))
                higgsLnKT.append(self.lnKT(genpart[higgs], t, delta, "tau2"))
                higgsLnZ.append(np.float32(math.log(z)))
                higgsLnDelta.append(np.float32(math.log(delta)))
                higgsLnKappa.append(np.float32(math.log(z * delta)))
                higgsLnPsi.append(self.Psi(genpart[higgs], t, "tau2"))
            else:
                tau2LnM.append(-99999.)
                tau2LnKT.append(-99999.)
                tau2LnZ.append(-99999.)
                tau2LnDelta.append(-99999.)
                tau2LnKappa.append(-99999.)
                tau2LnPsi.append(-99999.)
                higgsLnM.append(-99999.)
                higgsLnKT.append(-99999.)
                higgsLnZ.append(-99999.)
                higgsLnDelta.append(-99999.)
                higgsLnKappa.append(-99999.)
                higgsLnPsi.append(-99999.)

        self.out.fillBranch("SVFit_LundTau1M",          tau1LnM)
        self.out.fillBranch("SVFit_LundTau1KT",         tau1LnKT)
        self.out.fillBranch("SVFit_LundTau1Z",          tau1LnZ)
        self.out.fillBranch("SVFit_LundTau1Delta",      tau1LnDelta)
        self.out.fillBranch("SVFit_LundTau1Kappa",      tau1LnKappa)
        self.out.fillBranch("SVFit_LundTau1Psi",        tau1LnPsi)

        self.out.fillBranch("SVFit_LundTau2M",          tau2LnM)
        self.out.fillBranch("SVFit_LundTau2KT",         tau2LnKT)
        self.out.fillBranch("SVFit_LundTau2Z",          tau2LnZ)
        self.out.fillBranch("SVFit_LundTau2Delta",      tau2LnDelta)
        self.out.fillBranch("SVFit_LundTau2Kappa",      tau2LnKappa)
        self.out.fillBranch("SVFit_LundTau2Psi",        tau2LnPsi)
            
        self.out.fillBranch("SVFit_LundHiggsM",         higgsLnM)
        self.out.fillBranch("SVFit_LundHiggsKT",        higgsLnKT)
        self.out.fillBranch("SVFit_LundHiggsZ",         higgsLnZ)
        self.out.fillBranch("SVFit_LundHiggsDelta",     higgsLnDelta)
        self.out.fillBranch("SVFit_LundHiggsKappa",     higgsLnKappa)
        self.out.fillBranch("SVFit_LundHiggsPsi",       higgsLnPsi)
            

        return True
        
        
 # define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed
