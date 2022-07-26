# -*- coding: utf-8 -*-
"""
Created on Thu Jul  7 00:10:23 2022

@author: praya
"""

import sys
import os
import shutil
import numpy as np 
from subprocess import check_output
import matplotlib.pyplot as plt
import flopy.utils.binaryfile as bf
from numpy import absolute
from pandas import read_csv
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import RepeatedKFold
from xgboost import XGBRegressor
from xgboost import plot_tree
from xgboost import plot_importance
import pandas as pd
from matplotlib.pylab import rcParams
#rcParams['figure.figsize'] = 500,10
dataframe1 = pd.read_excel('Models5.xlsx')
data=dataframe1.values
d1=data[:,3:10]
d2=data[:,2]
X=d1
y=d2
model = XGBRegressor()
cv = RepeatedKFold(n_splits=10, n_repeats=3, random_state=10)
# evaluate model
scores = cross_val_score(model, X, y, scoring='r2', cv=cv, n_jobs=-1)
# force scores to be positive
scores = absolute(scores)
print('Mean MAE: %.3f (%.3f)' % (scores.mean(), scores.std()) )
model.fit(X,y)
def plot_tree(xgb_model, filename, rankdir='UT'):
    """
    Plot the tree in high resolution
    :param xgb_model: xgboost trained model
    :param filename: the pdf file where this is saved
    :param rankdir: direction of the tree: default Top-Down (UT), accepts:'LR' for left-to-right tree
    :return:
    """
    import xgboost as xgb
    import os
    gvz = xgb.to_graphviz(xgb_model, num_trees=xgb_model.best_iteration, rankdir=rankdir)
    _, file_extension = os.path.splitext(filename)
    format = file_extension.strip('.').lower()
    data = gvz.pipe(format=format)
    full_filename = filename
    with open(full_filename, 'wb') as f:
        f.write(data)
plot_tree(model,'xgboost_test_treeLR.pdf','LR')
plt.show()
plot_importance(model)