#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 13 18:43:35 2025

@author: alfonsocabezonvizoso
"""

'''
This script takes a CP sequence and performs ML predictions based on the best
models of classification and regression tasks.
The model used is a MPNN that has been trained on Coarse Grained graphs based on
the MA(R/S)TINI3 model.
'''


# =============================================================================
# Libraries
# =============================================================================

import itertools
import threading
import time
import sys
import torch
import dgl
from GraphGenerator import CP_Graph

# =============================================================================
# SCRIPT ARGUMENTS
# =============================================================================
import argparse
# Program description
parser = argparse.ArgumentParser(description =
    '''Cyclopeptide MD predictor \n
    This code takes a CP sequence and returns the prediction of the MD''' )

parser.add_argument( "-s", "--sequence", default = 'A',             #### NEEDED
    help = """Sequence of CP. Do not add D-amino acids""")
parser.add_argument( "-ff", "--forcefield", type = str, default = 'marstini3',
    help = """Forcefield.""")

args = parser.parse_args()


# =============================================================================
# FUNCTIONS
# =============================================================================
start = time.time()
done = False
def animate():
    '''
    Creates an animated loading

    Returns
    -------
    None.

    '''
    for c in itertools.cycle(['|', '/', '-', '\\']):
        if done:
            break
        sys.stdout.write('\rPredicting MD outcome ' + c)
        sys.stdout.flush()
        time.sleep(0.1)

t = threading.Thread(target=animate)
t.start()

# =============================================================================
# SCRIPT
# =============================================================================

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

# Create CP
CP = CP_Graph(args.sequence, args.forcefield)
g = CP.graph2dgl(pdb = f'{args.sequence}.pdb', itp = f'{args.sequence}.itp')
graph_batch = dgl.batch([g], ndata = ['h'], edata = ['length'])

# Load predictive models
model_reg = torch.load('model_reg_f1.pt', map_location=torch.device(device))
model_reg.eval() # set model to evaluation mode
model_class =  torch.load('model_Multiclass_f6.pt', map_location=torch.device(device))
model_class.eval() # set model to evaluation mode

# Class labels
map_class = {0 : "Lipophobic",
             1 : "Lipophilic",
             2 : "Intermediate"}

# =============================================================================
# RUN PREDICTIONS
# =============================================================================                                                             
                                                                                
with torch.no_grad():                                                       
    X = graph_batch.to(device)
    # Regression prediction                          
    pred_reg = model_reg(X) # Take reg predictions for test
    pred_reg =  torch.reshape(pred_reg, (-1,)) # reshape output
    pseudo_perm = pred_reg.numpy()[0]
    # Classification predictions
    logits_class = model_class(X) # Compute logits of the classes
    _, pred_class = torch.max(logits_class, 1) # Get predicted most likely cluster 
    class_num = pred_class.numpy()[0] # Get number of the cluster
    class_name = map_class[class_num]

# =============================================================================
# REPORT
# =============================================================================
done = True
form_seq = CP.format_seq(args.sequence)
print(f'\nReport:\nCP with sequence {form_seq} belongs to the {class_name} cluster with a predicted pseudo-permeability of {pseudo_perm:.2f}')
end = time.time()

print(f'--- TOTAL TIME: {end-start:.2f} seconds ---')