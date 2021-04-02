import os, sys
from os import path
import io
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import math
import numpy as np
import gzip
import json
from array import array
from importlib import import_module

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.tools import deltaPhi, deltaR, closest
from PhysicsTools.NanoAODTools.postprocessing.framework.treeReaderArrayTools import *
from rootpy.tree import Tree, TreeModel, IntCol, FloatArrayCol

class HiggsJSONProducer(Module):
    def __init__(self, processName):
        self.metBranchName = "MET"
        self.filename=processName + '.json.gz'

    def beginJob(self):
        pass
    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.fout = gzip.open(self.filename, 'w')

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.fout.close()

    def makeTLorentzVector(self, svfit, idx, isTau=""):
        v = ROOT.TLorentzVector()
        if isTau=="tau1":   v.SetPtEtaPhiM(svfit[idx].tau1Pt, svfit[idx].tau1Eta, svfit[idx].tau1Phi, svfit[idx].tau1Mass)
        elif isTau=="tau2": v.SetPtEtaPhiM(svfit[idx].tau2Pt, svfit[idx].tau2Eta, svfit[idx].tau2Phi, svfit[idx].tau2Mass)
        else:               v.SetPtEtaPhiM(svfit[idx].Pt, svfit[idx].Eta, svfit[idx].Phi, svfit[idx].Mass)
        return v

    def higg2json(self, svfit, idx):
        h = self.makeTLorentzVector(svfit, idx)
        t1 = self.makeTLorentzVector(svfit, idx, 'tau1')
        t2 = self.makeTLorentzVector(svfit, idx, 'tau2')

        
        j = [{'E':float(h.E()), 'px':float(h.Px()), 'py':float(h.Py()), 'pz':float(h.Pz())}]
        j.append({'E':float(t1.E()), 'px':float(t1.Px()), 'py':float(t1.Py()), 'pz':float(t1.Pz())})
        j.append({'E':float(t2.E()), 'px':float(t2.Px()), 'py':float(t2.Py()), 'pz':float(t2.Pz())})
        return j

    def analyze(self, event):
        ## Getting objects
        svfit     = Collection(event, "SVFit")
        svfitmet  = Collection(event, "SVFitMET")

        for idx in xrange(len(svfit)):
            if svfit[svfit[idx].Index].Pt > 0.:
                j = self.higg2json(svfit, svfit[idx].Index)
                self.fout.write((json.dumps(j)+'\n').encode('utf-8'))

        return True
        
        
 # define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed
