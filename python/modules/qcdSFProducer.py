#!/usr/bin/env python
import os, sys
import ROOT
import math
ROOT.PyConfig.IgnoreCommandLineOptions = True
from importlib import import_module
from os import system, environ
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.tools import deltaPhi, deltaR, closest
#import "$CMSSW_BASE/src/AnalysisTools/QuickRefold/interface/TObjectContainer.h"

class qcdSFProducer(Module): 
    def __init__(self, era):
	self.writeHistFile=True
	self.metBranchName="MET"
	self.debug = True
	if era == "2016":
        	self.correctionfile=os.environ["CMSSW_BASE"] + "/src/PhysicsTools/NanoSUSYTools/data/qcdJetRes/2016/qcdJetRespTailCorr.root"
	elif era == "2017":
        	self.correctionfile=os.environ["CMSSW_BASE"] + "/src/PhysicsTools/NanoSUSYTools/data/qcdJetRes/2017/qcdJetRespTailCorr.root"
	elif era == "2018":
        	self.correctionfile=os.environ["CMSSW_BASE"] + "/src/PhysicsTools/NanoSUSYTools/data/qcdJetRes/2018/qcdJetRespTailCorr.root"
        self.correctionhist="RespTailCorr"
      
    def beginJob(self,histFile=None,histDirName=None):
   	pass

    def endJob(self):
	pass 

    def loadhisto(self,filename,hname):
        file =ROOT.TFile.Open(filename)
        hist_ = file.Get(hname)
        hist_.SetDirectory(0)
        return hist_

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("qcdRespTailWeight","F")
        self.corrhist=self.loadhisto(self.correctionfile,self.correctionhist)

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def getQCDRespTailCorrector(self, jets, genJets, met):
      	MM = -1
      	ind = -1
      	flv = -1
      	resp = -1
	mmout = []
	for iG in range(0,len(genJets)):
		if iG > 4: break
        	if genJets[iG].pt == 0: break
        	fpt = -1
		for rJ in xrange(len(jets)):
			if not jets[rJ].Stop0l: continue
			if jets[rJ].genJetIdx == iG:
				fpt = jets[rJ].pt
				break
        	if fpt < 0: fpt = 9.5 #for the cases with no reco jet due to being below thresh
        	if(MM < 0 or abs(fpt - genJets[iG].pt) > MM):
			ind = iG
			resp =  fpt/genJets[iG].pt
			MM = abs(fpt - genJets[iG].pt)
			flv = abs(genJets[iG].partonFlavour)
    
	if ind >= 0:
		mmResp = resp
		mmInd = ind
		mmFlv = flv
	else:
		mmResp = -1
		mmInd = -1
		mmFlv = -1

	#mmout = [mmInd, mmResp, mmFlv]
	return mmResp, mmInd, mmFlv


    def analyze(self, event):
	jets      = Collection(event, 	"Jet")
	genjets   = Collection(event, 	"GenJet")
	met       = Object(event,     	self.metBranchName)
        stop0l 	  = Object(event, 	"Stop0l")

	mmResp, mmInd, mmFlv = self.getQCDRespTailCorrector(jets, genjets, met) 
	B_ratio = 1 if mmFlv != 5 else 2
	#B_ratio = 1 if abs(genjets[mmInd].partonFlavour) != 5 else 2

	qcdresptailweight = self.corrhist.GetBinContent(self.corrhist.GetXaxis().FindBin(mmResp),self.corrhist.GetYaxis().FindBin(B_ratio))
	
	if self.debug:
		print "xbin: ", self.corrhist.GetXaxis().FindBin(mmResp)
		print "ybin: ", self.corrhist.GetYaxis().FindBin(B_ratio)
		print "B ratio: ", B_ratio
		print "Response: ", mmResp
		print 'qcdrespTailWeight: ', qcdresptailweight 
        self.out.fillBranch("qcdRespTailWeight", qcdresptailweight)


        return True
