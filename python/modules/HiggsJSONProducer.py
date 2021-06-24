import os, sys
from os import path
import io
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import math
import numba
import numpy as np
import gzip
import json
from array import array
from importlib import import_module

from PhysicsTools.NanoSUSYTools.modules.HelperFunctions import *
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

def deltaRMatch(HiggsCandPdgId, HiggsCandEta, HiggsCandPhi, genTauDaughters_pdgId, genTauDaughters_eta, genTauDaughters_phi):

    matches = np.zeros(len(HiggsCandEta), dtype=int)
    matchesPdgId = np.zeros(len(HiggsCandEta), dtype=int)

    if len(genTauDaughters_eta):

        higgsEtaVals = np.array(np.meshgrid(HiggsCandEta, genTauDaughters_eta)).T.reshape(-1,2)
        higgsPhiVals = np.array(np.meshgrid(HiggsCandPhi, genTauDaughters_phi)).T.reshape(-1,2)
        higgsPdgIdVals = np.array(np.meshgrid(HiggsCandPdgId, genTauDaughters_pdgId)).T.reshape(-1,2)

        ## Using ufunc for vector operation
        deta = np.power(higgsEtaVals[:,0] - higgsEtaVals[:,1], 2)
        dPhi = higgsPhiVals[:,0] - higgsPhiVals[:,1]
        dR = np.sqrt((( abs(abs(dPhi)-np.pi)-np.pi )**2+(deta)**2)).reshape([-1,len(genTauDaughters_eta)/len(genTauDaughters_eta), len(genTauDaughters_eta)])
        ## Check PDG agreement
        matchPdgID = (higgsPdgIdVals[:,0] == abs(higgsPdgIdVals[:,1])).reshape([-1,len(genTauDaughters_eta)/len(genTauDaughters_eta), len(genTauDaughters_eta)])
        
        matches[dR.max(axis=2).min(axis=1) < 0.8] = 1
        matchesPdgId[matchPdgID.max(axis=2).min(axis=1)] = 1

    return matches, matchesPdgId

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

class HiggsJSONProducer(Module):
    def __init__(self, processName, match):
        self.metBranchName = "MET"
        self.processName = processName
        self.match = match
        self.filename='GenTau.json.gz'
        self.debug = False

    def beginJob(self):
        pass
    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.fout = gzip.open(self.filename, 'a')
        self.fgenout = gzip.open("genHiggs_"+self.filename, 'a')

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.fout.close()
        self.fgenout.close()

    def makeTLorentzVector(self, svfit, isTau=""):
        v = ROOT.TLorentzVector()
        if isTau=="tau1":   v.SetPtEtaPhiM(svfit.tau1Pt, svfit.tau1Eta, svfit.tau1Phi, svfit.tau1Mass)
        elif isTau=="tau2": v.SetPtEtaPhiM(svfit.tau2Pt, svfit.tau2Eta, svfit.tau2Phi, svfit.tau2Mass)
        elif isTau=="gen":  v.SetPtEtaPhiM(svfit.pt, svfit.eta, svfit.phi, svfit.mass)
        else:               v.SetPtEtaPhiM(svfit.Pt, svfit.Eta, svfit.Phi, svfit.Mass)
        return v

    def checkJSON(self, jin, jcheck):
        for j1, j in enumerate(jin):
            if self.debug: print("{0} and {1}".format(j, jcheck[0]))
            if j['E'] == jcheck[0]['E'] and j['px'] == jcheck[0]['px'] and j['py'] == jcheck[0]['py'] and j['pz'] == jcheck[0]['pz']:
                return False
        return True


    def higg2json(self, svfit, tau1GenMatch, tau2GenMatch, tau1GenMatchPdgId, tau2GenMatchPdgId, svfitmet, genHiggs = None):
        j = []
        isFirst = True
        for idx, sv in enumerate(svfit):
            if not sv.PassBaseline or not sv.PassLepton: continue
            if (tau1GenMatch[idx] and tau2GenMatch[idx]) or (tau1GenMatchPdgId[idx] and tau2GenMatchPdgId[idx]) or genHiggs == None:
                if isFirst: h = self.makeTLorentzVector(sv if genHiggs == None else genHiggs, "gen")
                t1 = self.makeTLorentzVector(sv, 'tau1')
                t2 = self.makeTLorentzVector(sv, 'tau2')
                
                if isFirst:
                    jTot = [{'E':float(h.E()), 'px':float(h.Px()), 'py':float(h.Py()), 'pz':float(h.Pz())}]
                    j.append(jTot[0])
                    nu1 = [{'E':float(sv.tau1nuE), 'px':float(sv.tau1nuPx), 'py':float(sv.tau1nuPy), 'pz':float(sv.tau1nuPz)}]
                    nu2 = [{'E':float(sv.tau2nuE), 'px':float(sv.tau2nuPx), 'py':float(sv.tau2nuPy), 'pz':float(sv.tau2nuPz)}]
                    j.append(nu1[0])
                    j.append(nu2[0])
                cand1 = [{'E':float(t1.E()), 'px':float(t1.Px()), 'py':float(t1.Py()), 'pz':float(t1.Pz())}]
                cand2 = [{'E':float(t2.E()), 'px':float(t2.Px()), 'py':float(t2.Py()), 'pz':float(t2.Pz())}]

                if self.debug: 
                    print("DM1: {0}, DM2: {1}".format(sv.tau1DM, sv.tau2DM))
                    print("cand1: {0}, cand2: {1}".format(cand1, cand2))
                fill1 = self.checkJSON(j, cand1)
                fill2 = self.checkJSON(j, cand2)
                if sv.tau1DM == -1. or sv.tau2DM == -1.:
                    fill1 = False
                    fill2 = False
                
                if fill1: j.append(cand1[0])
                if fill2: j.append(cand2[0])
                isFirst = False
        if self.debug: print("Output json: {0}".format(j))

        return j

    def tau2json(self, svfit, tau1GenMatch, tau1GenMatchPdgId, svfitmet, whichTau="tau1", genHiggs = None):
        j = []
        isFirst = True
        for idx, sv in enumerate(svfit):
            if not sv.PassBaseline or not sv.PassLepton: continue
            if (tau1GenMatch[idx]) or (tau1GenMatchPdgId[idx]) or genHiggs == None:
                if isFirst: h = self.makeTLorentzVector(sv if genHiggs == None else genHiggs, "gen")
                t1 = self.makeTLorentzVector(sv, whichTau)
                
                if isFirst:
                    jTot = [{'E':float(h.E()), 'px':float(h.Px()), 'py':float(h.Py()), 'pz':float(h.Pz())}]
                    j.append(jTot[0])
                    if whichTau == "tau1": nu1 = [{'E':float(sv.tau1nuE), 'px':float(sv.tau1nuPx), 'py':float(sv.tau1nuPy), 'pz':float(sv.tau1nuPz)}]
                    elif whichTau == "tau2": nu1 = [{'E':float(sv.tau2nuE), 'px':float(sv.tau2nuPx), 'py':float(sv.tau2nuPy), 'pz':float(sv.tau2nuPz)}]
                    j.append(nu1[0])
                cand1 = [{'E':float(t1.E()), 'px':float(t1.Px()), 'py':float(t1.Py()), 'pz':float(t1.Pz())}]

                if self.debug: 
                    print("DM1: {0}, cand1: {1}".format(sv.tau1DM, cand1))
                fill1 = self.checkJSON(j, cand1)
                if (sv.tau1DM == -1. and whichTau == "tau1") or (sv.tau2DM == -1. and whichTau == "tau2"):
                    fill1 = False
                
                if fill1: j.append(cand1[0])
                isFirst = False
        if self.debug: print("Output json: {0}".format(j))

        return j

    class TTreeReaderArrayWrapper:
        def __init__(self, ttarray):
            self.ttarray = ttarray

        def __iter__(self):
            for i in xrange(len(self.ttarray)):
                yield self.ttarray[i]
            return

    def HiggsGenMatch(self, event, HiggsCandEta, HiggsCandPhi, HiggsCandPdgId, match = 'GenPart'):
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
   
        matches, matchPdgID = deltaRMatch(HiggsCandPdgId, HiggsCandEta, HiggsCandPhi, genTauDaughters_pdgid, genTauDaughters_eta, genTauDaughters_phi)

        return matches, matchPdgID, genTauDaughters

    def recursiveFindHiggs(self, startIdx, targetPdgId, genpart):
        if startIdx < 0:
            return -1
    
        mom = genpart[startIdx].genPartIdxMother
    
        if mom < 0:
            return -1
        elif genpart[startIdx].pdgId == targetPdgId:
            return mom
        else:
            return self.recursiveFindHiggs(mom, targetPdgId, genpart)

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
        tau1GenMatch, tau1GenMatchPdgId, whichTau1 = self.HiggsGenMatch(event, tau1Eta, tau1Phi, tau1pdgIdTranslate, self.match)

        tau2Pt = np.fromiter(self.TTreeReaderArrayWrapper(event.SVFit_tau2Pt), dtype=float)
        tau2Eta = np.fromiter(self.TTreeReaderArrayWrapper(event.SVFit_tau2Eta), dtype=float)
        tau2Phi = np.fromiter(self.TTreeReaderArrayWrapper(event.SVFit_tau2Phi), dtype=float)
        tau2pdgId = np.fromiter(self.TTreeReaderArrayWrapper(event.SVFit_tau2pdgId), dtype=float)
        tau2DM = np.fromiter(self.TTreeReaderArrayWrapper(event.SVFit_tau2DM), dtype=float)
        tau2pdgIdTranslate = [tauDecayPdgId[TauDecay[t]] for t in tau2pdgId]
        tau2GenMatch, tau2GenMatchPdgId, whichTau2 = self.HiggsGenMatch(event, tau2Eta, tau2Phi, tau2pdgIdTranslate, self.match)

        tauIdx1 = []
        tauIdx2 = []
        higgsIdx = -1
        if self.debug: 
            print("tau1 pt: {0} eta: {1} pdgId: {2} DM: {3}".format(tau1Pt, tau1Eta, tau1pdgId, tau1DM))
            print("tau2 pt: {0} eta: {1} pdgId: {2} DM: {3}".format(tau2Pt, tau2Eta, tau2pdgId, tau2DM))
            print("tau1: {0} {1} {2} {3}".format(tau1GenMatch, whichTau1, tau1pdgIdTranslate, tau1GenMatchPdgId))
            print("tau2: {0} {1} {2} {3}".format(tau2GenMatch, whichTau2, tau2pdgIdTranslate, tau2GenMatchPdgId))
        for gD1, gD2 in zip(whichTau1, whichTau2):
            idx1 = genpart[gD1].genPartIdxMother
            idx2 = genpart[gD2].genPartIdxMother
            if idx1 >= 0: tauIdx1.append(idx1)
            if idx2 >= 0: tauIdx2.append(idx2)
            if idx1 >= 0:
                higgsIdx = genpart[idx1].genPartIdxMother
                if self.debug: print("Matched Higgs Index: {0} --> tau1: {1} pdgId: {3} tau2: {2} pdgId: {4}".format(genpart[idx1].genPartIdxMother, idx1, idx2, genpart[idx1].pdgId, genpart[idx2].pdgId))
                break
            elif idx2 >= 0:
                higgsIdx = genpart[idx2].genPartIdxMother
                if self.debug: print("Matched Higgs Index: {0} --> tau1: {1} pdgId: {3} tau2: {2} pdgId: {4}".format(genpart[idx1].genPartIdxMother, idx1, idx2, genpart[idx1].pdgId, genpart[idx2].pdgId))
                break

        j = self.tau2json(svfit, tau1GenMatch, tau1GenMatchPdgId, svfitmet, "tau1", genpart[idx1])
        if len(j) > 4: self.fout.write((json.dumps(j, sort_keys=False)+'\n').encode('utf-8'))
        j = self.tau2json(svfit, tau2GenMatch, tau2GenMatchPdgId, svfitmet, "tau2", genpart[idx2])
        if len(j) > 4: self.fout.write((json.dumps(j, sort_keys=False)+'\n').encode('utf-8'))
        if higgsIdx > 0: j = self.higg2json(svfit, tau1GenMatch, tau2GenMatch, tau1GenMatchPdgId, tau2GenMatchPdgId, svfitmet, genpart[higgsIdx])
        if len(j) > 4: self.fgenout.write((json.dumps(j, sort_keys=False)+'\n').encode('utf-8'))

        return True
        
        
 # define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed
