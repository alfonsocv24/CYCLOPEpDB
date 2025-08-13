#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 11 12:24:07 2024

@author: alfonsocabezonvizoso
"""

'''This script contains the GraphConv class, a custom graph convolutional
layer that will use node and edge features for message generation,
MLP for message passing, GRU for node update and Set2Set for readout
following the MPNN description by Gilmer et al as introduced in 
Neural Message Passing for Quantum Chemistry https://arxiv.org/abs/1704.01212'''


import torch.nn as nn
import torch.nn.functional as F
from dgl.nn.pytorch import NNConv, Set2Set

class MPNN(nn.Module):
    ''' This class performs message passing as the MPNN model and returns
    the updated node representation
    '''
    def __init__(self, node_in_feats : int, edge_in_feats : int = 1,
                 node_out_feats : int = 64, edge_hidden_feats : int = 128,
                 num_step_message_passing : int = 3):
        '''
        CPNN initializer. We will define the layers of our NN.
        
        Steps:
            1. Project node features linearly to reduce vector dimension.
            2. Define the edge metwork for message passing.
            3. Define the graph convolution using the edge network function.
            4. Define GRU as the update function

        Parameters
        ----------
        node_in_feats : int
            Size for the input node feature vector. 
        edge_in_feats : int, optional
            Size for the input edge feature vector. The default is 1 since we only
            have bond length.
        node_out_feats : int, optional
            Size for the updated node features. The default is 64.
        edge_hidden_feats : int, optional
            Size of edge hidden representation. The default is 128.
        num_step_message_passing : int, optional
            Number of message passing step (How far in the neighbours we are going).
            The default is 3.

        Returns
        -------
        None.

        '''
        super(MPNN, self).__init__()
        print(node_in_feats, node_out_feats)
        # Reduce node feature dimension
        self.project_node_feats = nn.Sequential(
            nn.Linear(node_in_feats, node_out_feats), nn.ReLU())
        
        self.num_step_message_passing = num_step_message_passing # number of message passing process
        # Define edge network to compute message
        edge_network = nn.Sequential(nn.Linear(edge_in_feats, edge_hidden_feats),
                                     nn.ReLU(), 
                                     nn.Linear(edge_hidden_feats, node_out_feats * node_out_feats))
        
        # Define graph convolution
        self.gnn = NNConv(in_feats = node_out_feats, out_feats= node_out_feats, 
                          edge_func = edge_network, aggregator_type='sum')
        
        # Define update function
        self.update = nn.GRU(node_out_feats, node_out_feats)
        
        
        
    def reset_parameters(self):
        '''Reinitializes model parameters'''
        self.project_node_feats[0].reset_parameters()
        self.gnn.reset_parameters()
        for layer in self.gnn.edge_func:
            if isinstance(layer, nn.Linear):
                layer.reset_parameters()
        self.update.reset_parameters()
    
    def forward(self, g, h, e):
        '''
        Function to launch the computation. Computes messages, performs message
        passing and updates node representations

        Parameters
        ----------
        g : DGLGraph
            DGLGraph for a batch of graphs.
        h : tensor (nodes_in_batch, node_in_feats)
            Input node feature tensor for the batch of graphs.
        e : tensor (edges_in_batch, edges_in_feats)
            Input edge feature tensor for the batch of graphs.

        Returns
        -------
        h : tensor of shape (nodes_in_batch, node_out_feat).
            Updated node feature tensor.

        '''
        h = self.project_node_feats(h) # Reduce node feature dimension
        ht = h.unsqueeze(0) # Adds dimension at specified position. Hidden state
        
        for _ in range(self.num_step_message_passing):
            h = F.relu(self.gnn(g, h, e)) # Apply convolution
            h, ht = self.update(h.unsqueeze(0), ht) # Updates and returns node features and final hidden state
            h = h.squeeze(0)
        return h
                

class GraphPredictor(nn.Module):
    ''' Graph Predictor based on a MPNN with Set2Set readout function and simple
    MLP for final regression/classification '''
    def __init__(self, node_in_feats : int, edge_in_feats : int = 1, 
                 node_out_feats : int = 64, edge_hidden_feats : int = 128,
                 ntasks : int = 1, nsteps_message_pass : int = 3,
                 nsteps_set2set : int = 3, nlayers_set2set : int = 2, 
                 task : str = 'Classification'):
        '''
        MPNN aaplied to classification and regression with set2set readout function
        Model introduced in Neural Message Passing for Quantum Chemistry.
        
        The initializer will define the layers of the model.
        
        Steps:
            1. Customized MPNN.
            2. Readout function (Set2Set)
            3. MLP for classification or regression
        

        Parameters
        ----------
        node_in_feats : int
            Size of the node feature vector.
        edge_in_feats : int, optional
            Size of the edge feature vector. Since we only have length, the default is 1.
        node_out_feats : int, optional
            Size for the output node feature vector. The default is 64.
        edge_hidden_feats : int, optional
            Size for the hidden edge representation. The default is 128.
        ntasks : int, optional
            Number of classes for classification, also output size. The default is 1.
        nsteps_message_pass : int, optional
            Number of message passing processes. The default is 3.
        nsteps_set2set : int, optional
            Number of set2set iterations. The default is 3.
        nlayers_set2set : int, optional
            Number of recurrent layers in set2set. The default is 2.
        task : str, optional
            Whether to perform Classification or Regression task. The default is Classification
        '''
        super(GraphPredictor, self).__init__()
        
        self.gnn = MPNN(node_in_feats = node_in_feats, node_out_feats = node_out_feats,
                        edge_in_feats = edge_in_feats, edge_hidden_feats = edge_hidden_feats,
                        num_step_message_passing = nsteps_message_pass) # MPNN framework defined previously
        # Use set2set as readout function
        self.readout = Set2Set(node_out_feats, nsteps_set2set, nlayers_set2set)
        
        # If we are doing classitication task, use MLP with sigmoid activation function
        # If we are doing regression, use MLP with no activation function in the last layer
        if task == 'Classification':
            if ntasks == 1:
                self.predictor = nn.Sequential(nn.Linear(2 * node_out_feats, node_out_feats),
                                                   nn.ReLU(),
                                                   nn.Dropout(0.2),
                                                   nn.Linear(node_out_feats, ntasks),
                                                   nn.Sigmoid())
            else:
                self.predictor = nn.Sequential(nn.Linear(2 * node_out_feats, node_out_feats),
                                                   nn.ReLU(),
                                                   nn.Dropout(0.2),
                                                   nn.Linear(node_out_feats, ntasks))
        if task == 'Regression':
            self.predictor = nn.Sequential(nn.Linear(2 * node_out_feats, node_out_feats),
                                               nn.ReLU(),
                                               nn.Linear(node_out_feats, ntasks))
        
    def forward(self, g):
        '''
        Performs Graph Classification/Regression

        Parameters
        ----------
        g : DGLGraph
            Batch of Graphs.
        h : torch.Tensor of shape (V, node_in_feats)
            Node feature tensor.
        e : torch.Tensor of shape (E, edge_in_feats)
            Edge feature tensor.

        Returns
        -------
        torch.Tensor of shape (G, ntasks).
            Prediction for the graphs in the batch. G is the number of graphs in batch

        '''
        h = g.ndata['h'].float()
        e = g.edata['length'].float()
        hT = self.gnn(g, h, e) # Perform Message passing and obtain final updated node features
        H = self.readout(g, hT) # Convert graph with updated node into single Tensor
        # H = self.dropout(H)
        return self.predictor(H)
            

class EarlyStopper:
    def __init__(self, patience : int = 1, min_delta : int = 0):
        self.patience = patience
        self.min_delta = min_delta
        self.counter = 0
        self.best_score = float('inf')
        self.best_accuracy = None
        self.best_model_state = None
        self.early_stop = False
        self.epoch = 0
        self.best_epoch = None

    def __call__(self, val_loss : float, val_acc : int, model):
        score = val_loss
        acc = val_acc
        self.epoch += 1
        if score < self.best_score:
            self.best_score = score
            self.best_accuracy = acc
            self.best_epoch = self.epoch
            self.counter = 0
            self.best_model_state = model.state_dict()
        elif score > (self.best_score + self.min_delta):
            self.counter += 1
            if self.counter >= self.patience:
                self.early_stop = True
        else:
            self.best_score = score
            self.best_accuracy = acc
            self.best_model_state = model.state_dict()
            self.best_epoch = self.epoch
            self.counter = 0
        
    def load_best_model(self, model):
        model.load_state_dict(self.best_model_state)
    
    def best_performance(self, task : str = 'Classification'):
        print(f'Best Validation loss: {self.best_score:.4f}. Achieved at epoch {self.best_epoch}')
        if task == 'Classification':
            print(f'Validation accuracy at epoch {self.best_epoch}: {self.best_accuracy * 100:.2f}%')
        elif task == 'Regression':
            print(f'Validation MAE at epoch {self.best_epoch}: {self.best_accuracy:.4f}')
        return self.best_epoch
    def reset(self):
        self.counter = 0
        self.best_score = float('inf')
        self.best_accuracy = None
        self.best_model_state = None
        self.early_stop = False
        self.epoch = 0
        self.best_epoch = None
        
        