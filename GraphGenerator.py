#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 15:55:42 2024

@author: alfonsocabezonvizoso
"""

'''This script contains the CYCLOPEp class, which contains all the functions involved
in the creation of the database and posterior processing of the data including 
graph generation of the CPs from the pdb'''


import MDAnalysis as mda
from typing import Tuple
import pandas as pd
pd.set_option('future.no_silent_downcasting', True)
import numpy as np
from biopandas.pdb import PandasPdb
import pickle
import networkx as nx
import matplotlib.pyplot as plt
import dgl
import os
os.environ["LOGURU_LEVEL"] = "ERROR"   # hides WARNING and below
from graphein.protein.visualisation import plotly_protein_structure_graph
from graphein.protein.graphs import label_node_id
from graphein.protein.graphs import initialise_graph_with_metadata
from graphein.protein.graphs import add_nodes_to_graph
from CP_PDB_Generator import Cyclopeptide


#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# SCRIPT ARGUMENTS
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
import argparse
# Program description
parser = argparse.ArgumentParser(description =
    '''Cyclopeptide Graph Generator \n
    Create the PDB and ITP files of a CP given its sequence.\n
    Then create its graph for prediction with ML''' )

parser.add_argument( "-s", "--sequence", default = 'A',             #### NEEDED
    help = """Sequence of CPs.""")
parser.add_argument( "-ff", "--forcefield", type = str, default = 'marstini3',
    help = """Forcefield.""")

args = parser.parse_args()


#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# FUNCTIONS
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# import graphein.molecule as gm


class CP_Graph(Cyclopeptide):
    def __init__(self, sequence, forcefield):
        # Initialize base class Cyclopeptide
        super().__init__(sequence, forcefield)
        # Assemble the CPs
        self.GenerateGeometry()
        
        # Generate the nanotube
        self.GenerateTopology()
        
        self.OptimizeStructure()
        self.sequence = self.format_seq(sequence)
        
    def format_seq(self, seq : str) -> str:
        '''
        This function takes the sequence of a CP and formats the string to add
        the d to the D-amino acids.

        Parameters
        ----------
        seq : str
            Sequence of the CP.

        Returns
        -------
        str
            Formatted sequence.

        '''
        aas = [aa for aa in seq]
        new_aas = [f'd{aa}'  if idx % 2 == 0 else aa for idx, aa in enumerate(aas)] # Converts even aa to d aa
        sequence = ''.join(new_aas)
        return sequence
        
    def bonds2pdb(self, pdb : str, itp : str):
        '''
        This function takes a pdb without conect section and the itp of the molecule
        and creates a pdb with the conect section

        Parameters
        ----------
        pdb : str.
            PDB file.
        itp : str
            ITP file of the topology of the molecule.

        Returns
        -------
        None.

        '''
        u = mda.Universe(itp, pdb, topology_format = 'ITP',
                         format = 'PDB') # Load itp and pdb of molecule to get bonds
        all_u = u.select_atoms('all')
        all_u.write(pdb, bonds = 'conect')
        
    def pdb2df(self, pdb : str, itp : str) -> Tuple['pd.core.frame.DataFrame', 'pd.core.frame.DataFrame', dict]:
        '''
        This function generates a pandas dataframe from a pdb file using biopandas
        Steps:
            1. Create universe of CP using pdb and itp
            2. Get bead type and their charges
            3. Create df from pdb
            4. Correct charges and add bead type to altloc column. Modify other columns
            5. Eliminate empty columns
            6. Create a node_id column using label_node_id function
            
        Parameters
        ----------
        pdb : str
            PDB file of CP.
        itp : str
            Topology .itp file of CP.

        Returns
        -------
        CP_df : 'pandas.core.frame.DataFrame'
            Formatted pandas df of the pdb.
        raw_CP_df : 'pandas.core.frame.DataFrame'
            Raw pandas df obtained from biopandas.
        mapping : dict.
            Dictionary of new labels for nodes based in beadtype:particle_number

        '''
        u = mda.Universe(itp , pdb, format = 'PDB' ,topology_format = 'ITP') # Create universe with molecule
        charges = u.atoms.charges # Get charge of each particle
        bead_type = u.atoms.types # Get bead type of each particle
        CP_df = PandasPdb().read_pdb(pdb) # Read pdb to pandas df
        CP_df = CP_df.df['ATOM'] # Within available dfs select ATOM (where our particles are)
        CP_df['charge'] = charges # Change charges to those of the itp
        CP_df['alt_loc'] = bead_type # Add bead type to alt_loc column
        CP_df['chain_id'] = 'CP' # Set Chain id name 
        CP_df['element_symbol'] = 'CG_bead'
        raw_CP_df = CP_df.copy() # Create a copy of the df before eliminate columns
        CP_df.replace("", float("NaN"), inplace = True) # Chanfe "" by Nan
        CP_df.dropna(axis=1, inplace = True) # eliminate NaN columns
        CP_df = label_node_id(CP_df, granularity = 'atom') # Adds a node_id column of form chain_id:residue_name:alt_loc:atom_name
        new_labels = [f'{bead}:{idx+1}' for idx, bead in enumerate(bead_type)]
        mapping = dict(zip(CP_df['node_id'].to_numpy(), new_labels))
        return CP_df, raw_CP_df, mapping
    
    def get_bond_dict(self, pdb : str, itp : str, length : int = 1.053) -> dict:
        '''
        This function takes a pdb and a itp of a CP and creates a dictionary for 
        the bonds storing the atoms involved and the distance
        
        Steps:
            1. Load molecule into Universe.
            2. Create a dict to store bonds in which keys are bonds numbered and
            the value is a list of [atom1, atom2, length] where atom1 and atom2 are atoms involved in bond
            3. Iterate over bonds to fill the dict

        Parameters
        ----------
        pdb : str
            Name of pdb file.
        itp : str
            Name of topology .itp file of the CP.
        length : int, optional
            Distance of the extra particles to BB. The default is 1.053.

        Returns
        -------
        dict
            dictionary with all bonds for the graph edges.

        '''
        u = mda.Universe(itp , pdb, 
                         format = 'PDB' ,topology_format = 'ITP') # Create universe with molecule
        bond_dict = {f'Bond{idx+1}' : [] for idx, bond in enumerate(u.bonds)} # Create dictionary to store bonds
        # Iterate over bonds and store the atom pairs of bonds and its length
        for idx, bond in enumerate(u.bonds):
            for atom in bond:
                bond_dict[f'Bond{idx+1}'].append(atom.ix)
            bond_dict[f'Bond{idx+1}'].append(bond.length())
            
        # Add bonds that correspond to extra particles
        for resid in u.residues:# Loop over residues
            atms_type = [] # Initialize atoms type container
            atms_idx = [] # Initialize atoms index container
            for atom in resid.atoms:# Loop over atoms in residue
                if 'CP' in atom.type:# It it is custom particle store them
                    atms_type.append(atom.type) # store atom type
                    atms_idx.append(atom.ix) # store index
# =============================================================================
#             Now we will get CPBB values and eliminate them from lists, this leads
#             to a more general code since it doesn't matter if CPCO or CPNH come
#             before CPBB, we will always have the CPBB-CPCO and CPBB-CPNH values'
# =============================================================================
            CPBB = atms_type.index('CPBB') # Get position of CPBB
            CPBB_idx = atms_idx[CPBB] # Get index of CPBB
            # Eliminate CPBB from both lists 
            atms_idx.pop(CPBB)
            atms_type.pop(CPBB)
            bond1 = f'CPBB:{atms_type[0]}:{resid.resid}' # Create key for dict based on resid
            bond2 = f'CPBB:{atms_type[1]}:{resid.resid}' # key for dict based on resid
            bond_dict[bond1] = [CPBB_idx, atms_idx[0], length] # Add value to dict
            bond_dict[bond2] = [CPBB_idx, atms_idx[1], length] # Add value to dict
        return bond_dict
    
    def gather_node_data(self, pdb : str, itp : str) -> dict:
        '''
        This function uses the information of the pdb and itp files of a CP as well as
        the information regarding the forcefield stored in bead_info.pkl and
        Interactions_marstini3.csv.
        STEPS:
            1. Load all the data that we need.
            2. Create a nested dictionary in which the keys are the nodes of the graph.
            The values will be dictionaries containing the info
            3. Loop over the atoms in the pdb and gather mass and charge info.
            4. Inside the loop, loop over the beads in the ff to gather info
            about epsilon and sigma values. This nested loop ensures always the same
            organization of the data
            5. Assign the info to the corresponding node

        Parameters
        ----------
        pdb : str
            PDB file of the CP.
        itp : str
            ITP topology file of CP.

        Returns
        -------
        features : dict
            DESCRIPTION.

        '''
        with open('bead_info.pkl', 'rb') as f:
            bead_info = pickle.load(f) # Load bead info of force field
        beads_in_ff = np.array(list(bead_info.keys()))
        cols_ints = ['Bead1', 'Bead2', 'Func', 'eps', 'sigm'] # Names for cols
        interactions = pd.read_csv('Interactions_marstini3.csv', index_col=False,
                            sep = '\s+', header = None, 
                            names = cols_ints) # Create Df for interactions
        # sys.exit()
        u = mda.Universe(itp , pdb, 
                         format = 'PDB' ,topology_format = 'ITP') # Create universe with molecule
        features = {f'{atom.type}:{atom.ix + 1}' : { 'h' :[] } for atom in u.atoms} # Create features dict with nodes names as keys
        # Loop over atoms to collect info
        for atom in u.atoms:
            coords = list(atom.position)
            charge = round(atom.charge, 1) # Get charge from itp
            mass = bead_info[atom.type]['Mass'] # Get mass from dict
            # Filter interaction df to get rows in which the current bead appears
            inters = interactions.loc[(interactions['Bead1'] == atom.type) | (interactions['Bead2'] == atom.type)]#, ['eps', 'sigm']]
            epsilons = [] # epsilon values container 
            sigmas = [] # sigma values container
# =============================================================================
#             Now we will go through the beads in the forcefield and find their
#             interactions with the current bead. By doing it in this way, we will
#             ensure that the interactions are all organized in the same way for 
#             each bead. The order is dictated by the order in which beads appear in 
#             the martini 3 ff
# =============================================================================
            for bead in beads_in_ff:
                if bead == atom.type:
                    # For self interations
                    filtered_inters = inters.loc[(inters['Bead1'] == bead) & (inters['Bead2'] == bead), ['eps', 'sigm']]
                    epsilons.append(filtered_inters['eps'].iloc[0])
                    sigmas.append(filtered_inters['sigm'].iloc[0])
                else:
                    filtered_inters = inters.loc[(inters['Bead1'] == bead) | (inters['Bead2'] == bead), ['eps', 'sigm']]
                    if len(filtered_inters.index) == 0:
                        # If there are no interactions, case of CPCO and CPNH with other particles, set a 0.0
                        epsilons.append(0.0)
                        sigmas.append(0.0)
                    else:
                        # Get values from file
                        epsilons.append(filtered_inters['eps'].iloc[0])
                        sigmas.append(filtered_inters['sigm'].iloc[0])
            # Add data to dict
            features[f'{atom.type}:{atom.ix + 1}']['h'].extend(coords)
            features[f'{atom.type}:{atom.ix + 1}']['h'].append(mass)
            features[f'{atom.type}:{atom.ix + 1}']['h'].append(charge)
            features[f'{atom.type}:{atom.ix + 1}']['h'].extend(epsilons)
            features[f'{atom.type}:{atom.ix + 1}']['h'].extend(sigmas)
        return features
    
    def clear_graph(self, graph : 'nx.classes.graph.Graph') -> 'nx.classes.graph.Graph':
        '''
        This function eliminates unnecessary attributes from CP graph. The attributes
        that we will eliminate are stored in attributes variable.
        
        Steps:
            1. Select attributes to remove.
            2. Loop over nodes and then over the created list of attribute
            3. Use the nested loop to remove attributes with del
            4. Return Graph

        Parameters
        ----------
        graph : nx.classes.graph.Graph
            Graph of a Cyclic Peptide.

        Returns
        -------
        graph : networkx.classes.graph.Graph
            Graph of the CP after eliminating unnecesary attributes.

        '''
        attributes = ['chain_id', 'residue_name', 'residue_number', 'atom_type',
                'element_symbol', 'b_factor', 'coords'] # Attributes we want to eliminate from node
        for node in list(graph.nodes()):
            for attribute in attributes:
                del graph.nodes[node][attribute]
        del graph.graph['pdb_code']
        return graph
    
    def feature2h(self, graph : 'nx.classes.graph.Graph') -> 'nx.classes.graph.Graph':
        '''
        Collapse all the node feratures into one called h.
        
        Steps:
            1. Get node list and attribute list.
            2. Nested loop to iterate over nodes and attributes.
            3. Append attributes to single list and remove original from node.
            4. Assign new list to feature in node

        Parameters
        ----------
        graph : 'networkx.classes.graph.Graph'
            Graph representing the CP.

        Returns
        -------
        graph : networkx.classes.graph.Graph
            Graph of the CP with the new node features.

        '''
        # Get node attribute names 
        with open('./scaler_features.pkl', 'rb') as f:
            scaler, _ = pickle.load(f)
        nodes = list(graph.nodes())
        for node in nodes:
            h = scaler.transform([graph.nodes[node]['h']])
            graph.nodes[node]['h'] = h[0]
        return graph
    
    def gather_features(self, pdbs : list, itps : list):
        '''
        Collapse all the node feratures into one called h.
        
        Steps:
            1. Get node list and attribute list.
            2. Nested loop to iterate over nodes and attributes.
            3. Append attributes to single list and remove original from node.
            4. Assign new list to feature in node

        Parameters
        ----------
        pdbs : list
            List of PDB files of the molecule
        itps : list
            List of ITP file of the molecule

        Returns
        -------
        None

        '''
        features = [] # Container for features
        # Iterate over list of pdbs
        counter = 0
        for idx, pdb in enumerate(pdbs):
            counter += 1
            itp = itps[idx]
            # CP_id, CP_name, sequence = self.infoDB(pdb) # Get info of CP from DataBase
            graph, sequence = self.pdb2graph(pdb, itp, ret = True)
            graph = graph.to_directed() # convert graph to directed
            # Get node attribute names 
            nodes = list(graph.nodes())
            for node in nodes:
                features.append(graph.nodes[node]['h'])
            print(f'Seqs remaining: {len(pdbs) - counter}')
        return features
            
    def pdb2graph(self, pdb : str, itp : str, ret : bool = False):
        '''
        This function takes a .pdb and .itp file of a CP and creates a graph from it.
        The nodes represent beads and have as features the mass, the charge and 
        the epsilon and sigma values of the interactions of the bead they represent
        in the marstini 3 forcefield.
        The edges are the bonds and the bond length is the only feature.
        The graphs are stored in separate .csv files for nodes and edges.
        
        Steps:
            1. Get CP info from CYCLOPEp database
            2. Create a pdb-based df with biopandas
            3. Initialize the graph using the df as metadata
            4. Add nodes.
            5. Add edges and its feature
            6. Rename nodes, compute features and add them.
            7. Convert nodes to df, apply format and save as ID_name_sequence_nodes.csv
            8. Convert edges to df and save as ID_name_sequence_edges.csv

        Parameters
        ----------
        pdb : str
            .pdb file of the CP.
        itp : str
            .itp of the CP.
        ret : bool, optional
            Whether to return the graph or not. The default is False.

        Returns
        -------
        CP : networkx.classes.graph.Graph
            Graph of the CP
        CP_id : int
            ID of the CP in the database
        CP_name : str
            name of the CP in the database
        sequence : str
            Sequence of the CP as stored in the database

        '''
        sequence = self.sequence
        # self.bonds2pdb(pdb , itp) # Add connect section to pdb
        CP_df, raw_CP_df, mapping = self.pdb2df(pdb, itp) # Creates a df from the pdb structure
        CP = initialise_graph_with_metadata(protein_df = CP_df, raw_pdb_df = raw_CP_df,
                                               path = pdb, granularity = 'atom') # Initialize graph from pdb
        CP = add_nodes_to_graph(CP) # Add nodes with the identifier
        nodes = list(CP.nodes()) # create node list
        bond_dict = self.get_bond_dict(pdb, itp) # Get dictionary with bond info
        # Iterate over bonds in dictionary
        for idx, bond in enumerate(bond_dict.keys()):
            features = bond_dict[bond] # Get features
            CP.add_edge(nodes[features[0]], nodes[features[1]],
                        length = [features[2]]) # Add edge to graph
        CP = nx.relabel_nodes(CP, mapping) # relabel nodes to easier ones
        CP = self.clear_graph(CP) # remove attributes we don't need
        features = self.gather_node_data(pdb, itp) # Create dict of dicts with new attributes for nodes
        nx.set_node_attributes(CP, features) # Add new attributes for nodes
        # for node in CP.nodes():
        #     CP.nodes[node]['coords'] = list(CP.nodes[node]['coords']) # Convert coordinates to list
        CP.graph['sequence_CP'] = sequence
        return CP, sequence
    
    def plot_graph(self, pdb : str, itp : str):
        '''
        This function creates a networkx graph from a pdb an itp and plots it.
        
        Steps:
            1. Create graph calling pdb2graph function.
            2. Modify nodes IDs so beads with same name have same color.
            3. Add a bond attribute so edges have all same color.
            4. Plot graph.

        Parameters
        ----------
        pdb : str
            Coordinate file of the CP.
        itp : str
            Topology file of the CP.

        Returns
        -------
        None.

        '''
        graph, _, _, _ = self.pdb2graph(pdb, itp, True) # Create graph
        nodes = list(graph.nodes()) # create node list
        for node in nodes:
            graph.nodes[node]['name'] = node.split(':')[0] # Adapts node name for visualization
        edges = list(graph.edges())
        for edge in edges:
            graph.edges[edge]['kind'] = 'bond' # Add bond attribute for representation
        print([ graph.nodes[node]['name'] for node in graph.nodes()])
        p = plotly_protein_structure_graph(graph,
                                           plot_title='Graph of CP', 
                                           figsize=(800, 800),
                                           colour_nodes_by='name',
                                           node_colour_map=plt.cm.tab20,
                                           colour_edges_by='kind',
                                           node_size_min=30)
        p.show()
        
    
    def graph2dgl(self, pdb : str, itp : str) -> 'dgl.heterograph.DGLGraph':
        '''
        This function uses a pdb and an itp of a CP to create a Graph in the
        Deep Graph Library.
        
        Steps:
            1. use infoDB function to obtain information of the CP.
            2. Use the pdb2graph function to create a graph of the CP in the
                networkx library.
            3. Convert graph to directed. Basically it doubles the number of
                edges. Necessary to transform the graph to dgl.
            4. Get node attributes and edge attributes. They are parameters
                for the conversion to dgl.
            5. Convert graph to dgl library using dgl.from_networkx.
            6. Save graph.

        Parameters
        ----------
        pdb : str
            Coordinate file of the CP.
        itp : str
            Topology file of the CP.

        Returns
        -------
        g : dgl graph.

        '''
        # from dgl.data.utils import save_graphs
        # CP_id, CP_name, sequence = self.infoDB(pdb) # Get info of CP from DataBase
        graph, sequence = self.pdb2graph(pdb, itp, ret = True)
        graph = graph.to_directed() # convert graph to directed
        # Get edge attribute name
        edges = list(graph.edges())
        edge_attrs = np.array(list(graph.edges[edges[0]].keys()))
        g = dgl.from_networkx(graph, node_attrs = ['h'], edge_attrs = edge_attrs)
        return g
        
        
