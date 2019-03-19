#!/usr/bin/env python 
import os, sys 
import ROOT 
ROOT.PyConfig.IgnoreCommandLineOptions = True 
from importlib import import_module 
import numpy as np 
import xgboost as xgb 
import logging 

class XGBHelper:     
    def __init__(self, model_file, var_list):         
        self.bst = xgb.Booster(params={'nthread': 1}, model_file=model_file)         
        self.var_list = var_list         
        logging.info('Load XGBoost model %s, input variables:\n  %s' % (model_file, str(var_list)))     

    def eval(self, inputs):         
        dmat = xgb.DMatrix(np.array([[inputs[k] for k in self.var_list]]), feature_names=self.var_list)         
        return self.bst.predict(dmat)[0]
