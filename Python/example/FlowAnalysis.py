# -*- coding: utf-8 -*-
"""
Created on Wed Jun 22 16:35:20 2016

@author: Alonyan
"""

from FlowCytometryTools import FCMeasurement
import os
import numpy as np
import pandas as pd

def flowstats(samples,ChanList):
        flowstats.means = samples[0][ChanList].mean()
        flowstats.medians = samples[0][ChanList].median()
        flowstats.stds = samples[0][ChanList].std()
        flowstats.counts = samples[0][ChanList].count()        
        flowstats.IDlist = [samples[0].ID[0:-4]]
        for i in np.arange(1,len(samples)):
            flowstats.means = pd.concat([flowstats.means, samples[i][ChanList].mean()], axis=1, ignore_index=True)
            flowstats.medians = pd.concat([flowstats.medians, samples[i][ChanList].median()], axis=1, ignore_index=True)
            flowstats.stds = pd.concat([flowstats.stds, samples[i][ChanList].std()], axis=1, ignore_index=True)
            flowstats.counts = pd.concat([flowstats.counts, samples[i][ChanList].count()], axis=1, ignore_index=True)
            flowstats.IDlist = flowstats.IDlist + [samples[i].ID[0:-4]]
        flowstats.means.columns = flowstats.IDlist
        flowstats.medians.columns = flowstats.IDlist
        flowstats.stds.columns = flowstats.IDlist
        flowstats.counts.columns = flowstats.IDlist
        return flowstats
        