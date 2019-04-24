#!/usr/bin/env python                                                                                                                                   
import os, sys
import ROOT
import math
import numpy as np
from array import array
ROOT.PyConfig.IgnoreCommandLineOptions = True
from importlib import import_module

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.tools import deltaPhi, deltaR, closest
from PhysicsTools.NanoAODTools.postprocessing.framework.treeReaderArrayTools import *
from rootpy.tree import Tree, TreeModel, IntCol, FloatArrayCol

class JetResSkim(Module):
    def __init__(self, era):
	self.era = era

    def beginJob(self):
        pass
    def endjob(self):
        pass
    def beginFile(self,inputFile,outputFile,inputTree,wrappedOutputTree):
	self.out = wrappedOutputTree
	self.out.branch("weight"	,"F")
	self.out.branch("genjetpt"	,"F")
	self.out.branch("genjeteta"	,"F")
	self.out.branch("recojetpt"	,"F")
	self.out.branch("genjetrank"	,"I")
	self.out.branch("flavor"	,"I")
	self.out.branch("rechempt"	,"F")
	self.out.branch("genhempt"	,"F")

    def PassEventFilter(self, flags):
        # https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2#2016_data
        passEventFilter = None

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 2016 ~~~~~
        if self.era == "2016":
            ## Common filters
            passEventFilter = flags.goodVertices and flags.HBHENoiseFilter and \
                    flags.HBHENoiseIsoFilter and flags.EcalDeadCellTriggerPrimitiveFilter \
                    and flags.BadPFMuonFilter and flags.BadChargedCandidateFilter # Post 2016 ICHEP
                    # and flags.BadPFMuonSummer16Filter and # flags.BadChargedCandidateSummer16Filter
            ## Only data

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 2017 ~~~~~
        if self.era == "2017" or self.era == "2018":
            ## Common filters
            passEventFilter = flags.goodVertices and flags.HBHENoiseFilter and \
                    flags.HBHENoiseIsoFilter and flags.EcalDeadCellTriggerPrimitiveFilter \
                    and flags.BadPFMuonFilter and flags.BadChargedCandidateFilter \
                    and flags.ecalBadCalibFilter  ## Need to double check whether is outdated
            ## Only data

        return passEventFilter

    def PassJetID(self, jets):
        # https://twiki.cern.ch/twiki/bin/view/CMS/JetID#Recommendations_for_13_TeV_2017
        # For 2016, loose and tight ID is the same : https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2016
        # For 2017, only tight ID available: https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2017
        # Select jet pt > 30GeV, which is used in jet ID study:
        # https://indico.cern.ch/event/592313/contributions/2398387/attachments/1384663/2106630/JetID_JMARmeeting_7_12_2016.pdf
        jetIDs = [j.jetId & 0b010 for j in jets if j.pt > 30]
        return (0 not in jetIDs)

    def analyze(self, event):
	jets        = Collection(event, "Jet")
	genjets     = Collection(event, "GenJet")
	flags       = Object(event,     "Flag")
	weight      = event.genWeight
	eventNum    = event.event

	PassJetID  = self.PassJetID(jets)
	PassFilter = self.PassEventFilter(flags) and PassJetID

	if PassFilter and PassJetID:
		for gJ in xrange(len(genjets)):
			gJet = genjets[gJ]
			if gJet.pt < 20: continue  
			
			rJet = 0
			for iR in xrange(len(jets)) :
				if jets[iR].genJetIdx != gJ: continue
				rJet = jets[iR]
				break
			
			self.out.fillBranch("weight",	weight)
			self.out.fillBranch("genjetpt", 	gJet.pt)
			self.out.fillBranch("genjeteta", 	gJet.eta)
			self.out.fillBranch("recojetpt", 	rJet.pt if rJet != 0 else 9.5)
			self.out.fillBranch("genjetrank", 	min(gJ, 250))
			self.out.fillBranch("flavor",	gJet.partonFlavour)
		
			if(gJet.eta > -2.8 and gJet.eta < -1.6 and gJet.phi >-1.37 and gJet.phi < -1.07):
				self.out.fillBranch("genhempt", gJet.pt)
			else:
				self.out.fillBranch("genhempt", 0)
			if rJet != 0:
				if(rJet.eta > -2.8 and rJet.eta < -1.6 and rJet.phi >-1.37 and rJet.phi < -1.07): 
					self.out.fillBranch("rechempt",rJet.pt)
				else:
					self.out.fillBranch("rechempt", 0)
			else:
				self.out.fillBranch("rechempt", 0)

			self.out.fill()	
	return True  
