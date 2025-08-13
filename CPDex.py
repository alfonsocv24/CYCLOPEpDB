#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  4 11:17:29 2022

@author: Fabián Suárez-Lestón
fabian.suarez.leston@usc.es
MD.USE Innovations &
CiQUS - Center for Research in Biological Chemistry and Molecular Materials &
Applied Physics Department, Faculty of Physics
University of Santiago de Compostela (USC), Gz, Spain
"""

import numpy as np

'''
 
------------------------------------------------------------------------------

 CPDex
 
 Information regarding the amino acids for the building of CPs and nanotubes
 
------------------------------------------------------------------------------

    This file contains 3 dictionaries:
        
    - FFs: dict of ( dict of str: str or int ) 
        
    Data regarding the forcefield used for the construction of the structure 
    and topology. Each entry corresponds to a forcefield available in CPBuilder. 
    All AA/UA forcefields currently available are defined under the key "AA".
   
    For each key, a dictionary with the following 4 entries must be specified:
        > "where": str
            "defined" means that the FF is included in Gromacs and tools
            such as pdb2gmx can be used.
            "itp" means that it is defined in an external ITP file and the
            topology needs to be build by the code.
        > "CG": str
            "yes" or "no", for CG or AA/UA FFs.
        > "backbone": str
            A selection in the MDAnalysis syntax for the elements of the
            backbone.
        > "distance": float
            A distance between the reference of adjacent residues in the
            CP for the initial structure.
    
            
    - Natural: dict of ( dict of str: list[ str or int or float ] or array[ float ] )
        
    Data regarding natural amino acids. Each key is the standard 1-letter 
    identifier of the amino acid, in capitals.
     
    Each entry contains:
        > "name": str
            The 3-letter identifier of the amino acid
            
        > FFs.keys(): dict
            A entry for each one in the FFs dictionary, storing information
            about the residue for building the structure/topology of the AA
            in each FF.
            
            The entries in the dictionary depends on the type of FF. Fields 
            marked with an asterisk are mandatory; otherwise must be present
            if the value of "where" in the FFs dictionary is "defined":
                
            >> * "beads": list[ str ]
                List of labels of the particles in the peptide in the 
                structure file.
              
            >> "mass": list[ float ]  
                List of masses of the particles. If None, the default value 
                in the FF file is used. It must follow the same order as the 
                "beads" entry. Values in g/mol
                
            >> "FF name": list[ str ]
                List of the name of the bead in the forcefield. It must 
                follow the same order as the "beads" entry.
                
            >> * "coords": array[ float ]
                Reference coordinates for the amino acid. The backbone should
                be aligned in the positive direction of the X axis (N->C) and
                with the reference particle/point in the origin. It must 
                follow the same order as the "beads" entry. Values in Å.
                
            >> "charge": list[ float ]
                List of charges of the particles. It must follow the same 
                order as the "beads" entry. Values in units of e
                
            >> "BB bonds": NOT IN USE
            
            >> "SC bonds": array
                List of bonds between elements in the side chain. Each row
                defines a bond and contains 4 elements: the indices in the 
                "beads" list of the bonded element, the bond length and the
                force constant.
                
            >> "constaints": array
                List of constraints in the amino acid. Each row defines a 
                constraint and contais 3 elements: the indices in the "beads"
                list of the constrained elements and the value of the distance
                between them.
                
            >> "BB angles": NOT IN USE
                
            >> "SC angles": array
                List of angles within the elements of the side chain. Each row
                defines an angle and contains 5 elements: the indices in the 
                "beads" list of the elements forming the angle, the value of
                the angle and the force constant.
                
            >> "BB dihedral": NOT IN USE
            
            >> "SC dihedral": array
                List of dihedrals within the elements of the side chain. Each
                row defines a dihedral and contains 6 elements: the indices in
                the "beads" list of the elements forming the dihedral, the value
                of the dihedral and the force constant.
                
            >> * "BB ref": list[ int ]
                List of the indices in the "beads" list corresponding to the
                elements in the backbone.
                
            >> * "BB-SC ref": list[ int ]
                List of the indices in the "beads" list correspondind to the
                BB->SC vector.
                
            >> "BB-SC angle": float
                Value of the angle between the BB->SC vector and the vector
                defined by the (i-1)th and i-th elements of the BB.
                
            >> "vsites3": dict of int: list[ float ]
                Virtual sites in the amino acid. The key is the index in the 
                "beads" list of the virtual site and the value is a list with
                the parameters of the virtual site. CHECK
                
            >> "exclusion": list[ int ]
                Indices in the "beads" list of the particles to generate
                exclusions.
                
            >> "int exclusion": bool
                CHECK
                    
            >> "vsitesn":
                CHECK
                   

    - Artificial: dict of ( dict of str: list[ str or int or float ] or array[ float ] )
        
    Data regarding artificial amino acids. Not in use; it must have the same 
    structure as the Natural dictionary.
       
'''

FFs = { 
       "AA": 
           {'where' : 'defined',
            'CG' : 'no',
            'backbone': 'backbone',
            'distance': 1.2},
       "martini22": 
           { 'where' : 'itp',
             'CG' : 'yes',
             'backbone' : 'name BB',
             'distance': 3.8},
       "marstini22": 
           { 'where' : 'itp',
             'CG' : 'yes',
             'backbone' : 'name BB',
             'distance': 3.8},
        "martini3":
            { 'where' : 'itp',
              'CG' : 'yes',
              'backbone' : 'name BB',
              'distance': 3.8},
        "marstini3": 
            { 'where' : 'itp',
              'CG' : 'yes',
              'backbone' : 'name BB',
              'distance': 3.8},
       "martini22p": 
           { 'where' : 'itp',
             'CG' : 'yes',
             'backbone' : 'name BB',
             'distance': 3.8}, 
           }


Natural = {
    # Arginine
    "R" :
        { "name": "ARG",
          "AA": {
              "beads": [ "N", "H", "CA", "CB", "CG", "CD", "NE", "HE", "CZ", "NH1", "HH11", "HH12", "NH2", "HH21", "HH22", "C", "O" ],
              "coords": np.array([[ -7.070e-01, -1.225e+00,  5.000e-01 ],
                                  [ -7.070e-01, -1.225e+00,  1.500e+00 ],
                                  [  0.000e+00,  0.000e+00,  0.000e+00 ],
                                  [  1.520e+00, -0.000e+00,  5.070e-01 ],
                                  [  2.133e+00,  1.391e+00,  5.070e-01 ],
                                  [  3.653e+00,  1.337e+00,  4.870e-01 ],
                                  [  4.240e+00,  2.673e+00,  5.040e-01 ],
                                  [  3.615e+00,  3.454e+00,  5.480e-01 ],
                                  [  5.545e+00,  2.931e+00,  4.680e-01 ],
                                  [  6.429e+00,  1.938e+00,  4.110e-01 ],
                                  [  6.115e+00,  9.880e-01,  3.960e-01 ],
                                  [  7.408e+00,  2.142e+00,  3.840e-01 ],
                                  [  5.968e+00,  4.190e+00,  4.880e-01 ],
                                  [  5.307e+00,  4.939e+00,  5.320e-01 ],
                                  [  6.948e+00,  4.390e+00,  4.610e-01 ],
                                  [ -7.160e-01,  1.241e+00,  5.070e-01 ],
                                  [ -7.160e-01,  1.241e+00,  1.737e+00 ]]),
              "BB ref": [ 0, 2, 22, 23 ],
              "SC ref": [ 3, 4 ],
              },
          "martini22":  {
              "beads": [ "BB", "SC1", "SC2" ],
              "mass": [ None, None, None ],
              "FF name": [ "P5", "N0", "Qd" ],
              "coords": np.array([[ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00],
                                  [ 3.21984253e+00,  7.14133012e-08, -8.83578968e-08],
                                  [ 6.92256500e+00, -1.36323624e-01,  5.51608469e-02]]),
              "charge": [ 0, 0, 1 ],
              "BB bonds": [],
              "SC bonds": np.array([[ 0, 1, 0.33, 5000 ],
                                    [ 1, 2, 0.34, 5000 ]]),
              "constraints": [],
              "BB angles": [],
              "SC angles": np.array([[ 0, 1, 2, 180, 25 ]]),
              "BB dihedral": [],
              "SC dihedral": [],
              "BB ref": [ 0 ],
              "BB-SC ref": [ 0, 1 ],
              "BB-SC angle": 100,
              "vsites3": { },
              "exclusion": [ ],
              "int exclusion": False,
              "vsitesn": []
              },
          "marstini22": {
              "beads": [ "BB", "SC1", "SC2", "H", "O" ],
              "mass": [ None, None, None, 0.00, 0.00 ],
              "FF name": [ "CPBB", "N0", "Qd", "CPNH", "CPCO" ],
              "coords": np.array([[ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00],
                                  [ 3.21984253e+00,  7.14133012e-08, -8.83578968e-08],
                                  [ 6.92256500e+00, -1.36323624e-01,  5.51608469e-02],
                                  [ 0.00000000e+00,  3.50000000e-01,  1.05300000e+00],
                                  [ 0.00000000e+00, -3.50000000e-01,  1.05300000e+00]]),
              "charge": [ 0, 0, 1, 0.4, -0.4 ],
              "BB bonds": [],
              "SC bonds": np.array([[ 0, 1, 0.33, 5000 ],
                                    [ 1, 2, 0.34, 5000 ]]),
              "constraints": [],
              "BB angles": [],
              "SC angles": np.array([[ 0, 1, 2, 180, 25 ]]),
              "BB dihedral": [],
              "SC dihedral": [],
              "BB ref": [ 0 ],
              "BB-SC ref": [ 0, 1 ],
              "BB-SC angle": 100,
              "vsites3": { 3: [ 0.35, 0.00, 1.053 ],
                           4: [ 0.00, 0.35, 1.053 ] },
              "exclusion": [ 3, 4 ],
              "int exclusion": False,
              "vsitesn": []
              },
          "martini3": {
              "beads": [ "BB", "SC1", "SC2" ],
              "mass": [ None, None, None ],
              "FF name": [ "P2", "SC3", "SQ3p"],
              "coords": np.array([[ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00],
                                  [ 3.21984253e+00,  7.14133012e-08, -8.83578968e-08],
                                  [ 6.92256500e+00, -1.36323624e-01,  5.51608469e-02]]),
              "charge": [ 0, 0, 1 ],
              "BB bonds": [],
              "SC bonds": np.array([[ 0, 1, 0.33, 5000 ],
                                    [ 1, 2, 0.38, 5000 ]]),
              "constraints": [],
              "BB angles": [],
              "SC angles": np.array([[ 0, 1, 2, 180, 25 ]]),
              "BB dihedral": [],
              "SC dihedral": [],
              "BB ref": [ 0 ],
              "BB-SC ref": [ 0, 1 ],
              "BB-SC angle": 100,
              "vsites3": {},
              "exclusion": [ ],
              "int exclusion": False,
              "vsitesn": []
              },
          "marstini3": {
              "beads": [ "BB", "SC1", "SC2", "H", "O" ],
              "mass": [ None, None, None, 0.00, 0.00 ],
              "FF name": [ "CPBB", "SC3", "SQ3p", "CPNH", "CPCO" ],
              "coords": np.array([[ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00],
                                  [ 3.21984253e+00,  7.14133012e-08, -8.83578968e-08],
                                  [ 6.92256500e+00, -1.36323624e-01,  5.51608469e-02],
                                  [ 0.00000000e+00,  3.50000000e-01,  1.05300000e+00],
                                  [ 0.00000000e+00, -3.50000000e-01,  1.05300000e+00]]),
              "charge": [ 0, 0, 1, 0.4, -0.4 ],
              "BB bonds": [],
              "SC bonds": np.array([[ 0, 1, 0.33, 5000 ],
                                    [ 1, 2, 0.38, 5000 ]]),
              "constraints": [],
              "BB angles": [],
              "SC angles": np.array([[ 0, 1, 2, 180, 25 ]]),
              "BB dihedral": [],
              "SC dihedral": [],
              "BB ref": [ 0 ],
              "BB-SC ref": [ 0, 1 ],
              "BB-SC angle": 100,
              "vsites3": { 3: [ 0.35, 0.00, 1.053 ],
                           4: [ 0.00, 0.35, 1.053 ] },
              "exclusion": [ 3, 4 ],
              "int exclusion": False,
              "vsitesn": []
              },
          "martini22p": {
              
              } 
        },
        
    # Histidine
    "H" :
        { "name": "HIS",
          "AA" : {
              "beads": [ "N", "H", "CA", "CB", "CG", "ND1", "CD2", "HD2", "CE1", "HE1", "NE2", "HE2", "C", "O" ],
              "coords": np.array([[ -7.070e-01, -1.225e+00,  5.000e-01 ],
                                  [ -7.070e-01, -1.225e+00,  1.500e+00 ],
                                  [  0.000e+00,  0.000e+00,  0.000e+00 ],
                                  [  1.520e+00,  0.000e+00,  5.070e-01 ],
                                  [  2.120e+00,  1.364e+00,  5.070e-01 ],
                                  [  3.482e+00,  1.583e+00,  5.080e-01 ],
                                  [  1.535e+00,  2.581e+00,  5.070e-01 ],
                                  [  5.500e-01,  2.756e+00,  5.060e-01 ],
                                  [  3.696e+00,  2.885e+00,  5.080e-01 ],
                                  [  4.593e+00,  3.327e+00,  5.090e-01 ],
                                  [  2.503e+00,  3.521e+00,  5.070e-01 ],
                                  [  2.363e+00,  4.512e+00,  5.060e-01 ],
                                  [ -7.160e-01,  1.241e+00,  5.070e-01 ],
                                  [ -7.160e-01,  1.241e+00,  1.737e+00 ]]),
              "BB ref": [ 0, 2, 12, 13 ],
              "SC ref": [ 3 ],
              },
          "martini22":  {
              "beads": [ "BB", "SC1", "SC2", "SC3"],
              "mass": [ None, None, None, None ],
              "FF name": [ "P5", "SC4", "SP1", "SP1"],
              "coords": np.array([[ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00],
                                  [ 2.85862120e+00, -1.49837623e-07,  2.28756755e-07],
                                  [ 4.67208271e+00, -7.59461588e-01, -1.40423356e+00],
                                  [ 4.53404185e+00, -1.59507914e+00,  1.49562690e-01]]),
              "charge": [ 0, 0, 0, 0 ],
              "BB bonds": [],
              "SC bonds": np.array([[ 0, 1, 0.32, 7500 ]]),
              "constraints": np.array([[ 1, 2, 0.27 ],
                                       [ 1, 3, 0.27 ],
                                       [ 2, 3, 0.27 ]]),
              "BB angles": [],
              "SC angles": np.array([[ 0, 1, 2, 150, 50 ],
                                     [ 0, 1, 3, 150, 50 ]]),
              "BB dihedral": [],
              "SC dihedral": np.array([[ 0, 2, 3, 1, 0, 50 ]]),
              "BB ref": [ 0 ],
              "BB-SC ref": [ 0, 1 ],
              "BB-SC angle": 100,
              "vsites3": { },
              "exclusion": [ ],
              "int exclusion": False,
              "vsitesn": []
              },
          "marstini22": {
              "beads": [ "BB", "SC1", "SC2", "SC3", "H", "O"],
              "mass": [ None, None, None, None, 0.00, 0.00 ],
              "FF name": [ "CPBB", "SC4", "SP1", "SP1", "CPNH", "CPCO"],
              "coords": np.array([[ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00],
                                  [ 2.85862120e+00, -1.49837623e-07,  2.28756755e-07],
                                  [ 4.67208271e+00, -7.59461588e-01, -1.40423356e+00],
                                  [ 4.53404185e+00, -1.59507914e+00,  1.49562690e-01],
                                  [ 0.00000000e+00,  3.50000000e-01,  1.05300000e+00],
                                  [ 0.00000000e+00, -3.50000000e-01,  1.05300000e+00]]),
              "charge": [ 0, 0, 0, 0, 0.4, -0.4 ],
              "BB bonds": [],
              "SC bonds": np.array([[ 0, 1, 0.32, 7500 ]]),
              "constraints": np.array([[ 1, 2, 0.27 ],
                                       [ 1, 3, 0.27 ],
                                       [ 2, 3, 0.27 ]]),
              "BB angles": [],
              "SC angles": np.array([[ 0, 1, 2, 150, 50 ],
                                     [ 0, 1, 3, 150, 50 ]]),
              "BB dihedral": [],
              "SC dihedral": np.array([[ 0, 2, 3, 1, 0, 50 ]]),
              "BB ref": [ 0 ],
              "BB-SC ref": [ 0, 1 ],
              "BB-SC angle": 100,
              "vsites3": { 4: [ 0.35, 0.00, 1.053 ],
                            5: [ 0.00, 0.35, 1.053 ] },
              "exclusion": [ 4, 5 ],
              "int exclusion": False,
              "vsitesn": []
              },
          "martini3": {
              "beads": [ "BB", "SC1", "SC2", "SC3"],
              "mass": [ None, None, None, None ],
              "FF name": [ "P2", "TC4", "TN6d", "TN5a"],
              "coords": np.array([[ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00],
                                  [ 2.85862120e+00, -1.49837623e-07,  2.28756755e-07],
                                  [ 4.67208271e+00, -7.59461588e-01, -1.40423356e+00],
                                  [ 4.53404185e+00, -1.59507914e+00,  1.49562690e-01]]),
              "charge": [ 0, 0, 0, 0],
              "BB bonds": [],
              "SC bonds": np.array([[ 0, 1, 0.336, 7500 ]]),
              "constraints": np.array([[ 1, 2, 0.32 ],
                                       [ 1, 3, 0.30 ],
                                       [ 2, 3, 0.27 ]]),
              "BB angles": [],
              "SC angles": np.array([[ 0, 1, 2, 120, 50 ],
                                     [ 0, 1, 3, 120, 50 ]]),
              "BB dihedral": [],
              "SC dihedral": [],
              "BB ref": [ 0 ],
              "BB-SC ref": [ 0, 1 ],
              "BB-SC angle": 100,
              "vsites3": {},
              "exclusion": [ ],
              "int exclusion": True,
              "vsitesn": []
              },
          "marstini3": {
              "beads": [ "BB", "SC1", "SC2", "SC3", "H", "O"],
              "mass": [ None, None, None, None, 0.00, 0.00 ],
              "FF name": [ "CPBB", "TC4", "TN6d", "TN5a", "CPNH", "CPCO"],
              "coords": np.array([[ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00],
                                  [ 2.85862120e+00, -1.49837623e-07,  2.28756755e-07],
                                  [ 4.67208271e+00, -7.59461588e-01, -1.40423356e+00],
                                  [ 4.53404185e+00, -1.59507914e+00,  1.49562690e-01],
                                  [ 0.00000000e+00,  3.50000000e-01,  1.05300000e+00],
                                  [ 0.00000000e+00, -3.50000000e-01,  1.05300000e+00]]),
              "charge": [ 0, 0, 0, 0, 0.4, -0.4 ],
              "BB bonds": [],
              "SC bonds": np.array([[ 0, 1, 0.336, 7500 ]]),
              "constraints": np.array([[ 1, 2, 0.32 ],
                                       [ 1, 3, 0.30 ],
                                       [ 2, 3, 0.27 ]]),
              "BB angles": [],
              "SC angles": np.array([[ 0, 1, 2, 120, 50 ],
                                     [ 0, 1, 3, 120, 50 ]]),
              "BB dihedral": [],
              "SC dihedral": [],
              "BB ref": [ 0 ],
              "BB-SC ref": [ 0, 1 ],
              "BB-SC angle": 100,
              "vsites3": { 4: [ 0.35, 0.00, 1.053 ],
                            5: [ 0.00, 0.35, 1.053 ] },
              "exclusion": [ 4, 5 ],
              "int exclusion": True,
              "vsitesn": []
              },
          "martini22p": {
              
              } 
        },
        
    # Lysine
    "K" :
        { "name": "LYS",
          "AA" : {
              "beads": [ "N", "H", "CA", "CB", "CG", "CD", "CE", "NZ", "HZ1", "HZ2", "HZ3", "C", "O" ],
              "coords": np.array([[ -7.070e-01, -1.225e+00,  5.000e-01 ],
                                  [ -7.070e-01, -1.225e+00,  1.500e+00 ],
                                  [  0.000e+00,  0.000e+00,  0.000e+00 ],
                                  [  1.521e+00,  0.000e+00,  5.070e-01 ],
                                  [  2.134e+00,  1.390e+00,  5.070e-01 ],
                                  [  3.653e+00,  1.336e+00,  4.600e-01 ],
                                  [  4.241e+00,  2.672e+00,  4.690e-01 ],
                                  [  5.545e+00,  2.927e+00,  4.430e-01 ],
                                  [  5.692e+00,  3.916e+00,  4.560e-01 ],
                                  [  5.981e+00,  2.515e+00,  1.242e+00 ],
                                  [  5.944e+00,  2.542e+00, -3.900e-01 ],
                                  [ -7.160e-01,  1.241e+00,  5.070e-01 ],
                                  [ -7.160e-01,  1.241e+00,  1.737e+00 ]]),
              "BB ref": [ 0, 2, 11, 12 ],
              "SC ref": [ 3 ],
              },
          "martini22":  {
              "beads": [ "BB", "SC1", "SC2" ],
              "mass": [ None, None, None ],
              "FF name": [ "P5", "C3", "Qd" ],
              "coords": np.array([[ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00],
                                  [ 3.50092902e+00,  7.30257805e-08,  1.41067297e-07],
                                  [ 6.57851381e+00, -6.16725423e-01,  8.90774417e-02]]),
              "charge": [ 0, 0, 1 ],
              "BB bonds": [],
              "SC bonds": np.array([[ 0, 1, 0.33, 5000 ],
                                    [ 1, 2, 0.28, 5000 ]]),
              "constraints": [],
              "BB angles": [],
              "SC angles": np.array([[ 0, 1, 2, 180, 25 ]]),
              "BB dihedral": [],
              "SC dihedral": [],
              "BB ref": [ 0 ],
              "BB-SC ref": [ 0, 1 ],
              "BB-SC angle": 100,
              "vsites3": { },
              "exclusion": [ ],
              "int exclusion": False,
              "vsitesn": []
              },
          "marstini22": {
              "beads": [ "BB", "SC1", "SC2", "H", "O" ],
              "mass": [ None, None, None, 0.00, 0.00 ],
              "FF name": [ "CPBB", "C3", "Qd", "CPNH", "CPCO" ],
              "coords": np.array([[ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00],
                                  [ 3.50092902e+00,  7.30257805e-08,  1.41067297e-07],
                                  [ 6.57851381e+00, -6.16725423e-01,  8.90774417e-02],
                                  [ 0.00000000e+00,  3.50000000e-01,  1.05300000e+00],
                                  [ 0.00000000e+00, -3.50000000e-01,  1.05300000e+00]]),
              "charge": [ 0, 0, 1, 0.4, -0.4 ],
              "BB bonds": [],
              "SC bonds": np.array([[ 0, 1, 0.33, 5000 ],
                                    [ 1, 2, 0.28, 5000 ]]),
              "constraints": [],
              "BB angles": [],
              "SC angles": np.array([[ 0, 1, 2, 180, 25 ]]),
              "BB dihedral": [],
              "SC dihedral": [],
              "BB ref": [ 0 ],
              "BB-SC ref": [ 0, 1 ],
              "BB-SC angle": 100,
              "vsites3": { 3: [ 0.35, 0.00, 1.053 ],
                            4: [ 0.00, 0.35, 1.053 ] },
              "exclusion": [ 3, 4 ],
              "int exclusion": False,
              "vsitesn": []
              },
          "martini3": {
              "beads": [ "BB", "SC1", "SC2"],
              "mass": [ None, None, None ],
              "FF name": [ "P2", "SC3", "SQ4p"],
              "coords": np.array([[ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00],
                                  [ 3.50092902e+00,  7.30257805e-08,  1.41067297e-07],
                                  [ 6.57851381e+00, -6.16725423e-01,  8.90774417e-02]]),
              "charge": [ 0, 0, 1],
              "BB bonds": [],
              "SC bonds": np.array([[ 0, 1, 0.33, 5000 ],
                                    [ 1, 2, 0.36, 5000 ]]),
              "constraints": [],
              "BB angles": [],
              "SC angles": np.array([[ 0, 1, 2, 180, 25 ]]),
              "BB dihedral": [],
              "SC dihedral": [],
              "BB ref": [ 0 ],
              "BB-SC ref": [ 0, 1 ],
              "BB-SC angle": 100,
              "vsites3": { },
              "exclusion": [ ],
              "int exclusion": False,
              "vsitesn": []
              },
          "marstini3": {
              "beads": [ "BB", "SC1", "SC2", "H", "O" ],
              "mass": [ None, None, None, 0.00, 0.00 ],
              "FF name": [ "CPBB", "SC3", "SQ4p", "CPNH", "CPCO" ],
              "coords": np.array([[ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00],
                                  [ 3.50092902e+00,  7.30257805e-08,  1.41067297e-07],
                                  [ 6.57851381e+00, -6.16725423e-01,  8.90774417e-02],
                                  [ 0.00000000e+00,  3.50000000e-01,  1.05300000e+00],
                                  [ 0.00000000e+00, -3.50000000e-01,  1.05300000e+00]]),
              "charge": [ 0, 0, 1, 0.4, -0.4 ],
              "BB bonds": [],
              "SC bonds": np.array([[ 0, 1, 0.33, 5000 ],
                                    [ 1, 2, 0.36, 5000 ]]),
              "constraints": [],
              "BB angles": [],
              "SC angles": np.array([[ 0, 1, 2, 180, 25 ]]),
              "BB dihedral": [],
              "SC dihedral": [],
              "BB ref": [ 0 ],
              "BB-SC ref": [ 0, 1 ],
              "BB-SC angle": 100,
              "vsites3": { 3: [ 0.35, 0.00, 1.053 ],
                            4: [ 0.00, 0.35, 1.053 ] },
              "exclusion": [ 3, 4 ],
              "int exclusion": False,
              "vsitesn": []
              },
          "martini22p": {
              
              } 
        },
        
    # Aspartic acid
    "D" :
        { "name": "ASP",
          "AA" : {
              "beads": [ "N", "H", "CA", "CB", "CG", "OD1", "OD2", "C", "O" ],
              "coords": np.array([[ -7.070e-01, -1.225e+00,  5.000e-01 ],
                                  [ -7.070e-01, -1.225e+00,  1.500e+00 ],
                                  [  0.000e+00,  0.000e+00,  0.000e+00 ],
                                  [  1.520e+00, -0.000e+00,  5.070e-01 ],
                                  [  2.115e+00,  1.398e+00,  5.070e-01 ],
                                  [  3.358e+00,  1.534e+00,  5.070e-01 ],
                                  [  1.332e+00,  2.374e+00,  5.070e-01 ],
                                  [ -7.160e-01,  1.241e+00,  5.070e-01 ],
                                  [ -7.160e-01,  1.241e+00,  1.737e+00 ]]),
              "BB ref": [ 0, 2, 7, 8 ],
              "SC ref": [ 3 ],
              },
          "martini22":  {
              "beads": [ "BB", "SC1" ],
              "mass": [ None, None ],
              "FF name": [ "P5", "Qa" ],
              "coords": np.array([[ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00],
                                  [ 3.52791575e+00,  3.98701858e-08, -3.49091158e-08]]),
              "charge": [ 0, -1 ],
              "BB bonds": [],
              "SC bonds": np.array([[ 0, 1, 0.32, 7500 ]]),
              "constraints": [],
              "BB angles": [],
              "SC angles": [],
              "BB dihedral": [],
              "SC dihedral": [],
              "BB ref": [ 0 ],
              "BB-SC ref": [ 0, 1 ],
              "BB-SC angle": 100,
              "vsites3": { },
              "exclusion": [ ],
              "int exclusion": False,
              "vsitesn": []
              },
          "marstini22": {
              "beads": [ "BB", "SC1", "H", "O" ],
              "mass": [ None, None, 0.00, 0.00 ],
              "FF name": [ "CPBB", "Qa", "CPNH", "CPCO" ],
              "coords": np.array([[ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00],
                                  [ 3.52791575e+00,  3.98701858e-08, -3.49091158e-08],
                                  [ 0.00000000e+00,  3.50000000e-01,  1.05300000e+00],
                                  [ 0.00000000e+00, -3.50000000e-01,  1.05300000e+00]]),
              "charge": [ 0, -1, 0.4, -0.4 ],
              "BB bonds": [],
              "SC bonds": np.array([[ 0, 1, 0.32, 7500 ]]),
              "constraints": [],
              "BB angles": [],
              "SC angles": [],
              "BB dihedral": [],
              "SC dihedral": [],
              "BB ref": [ 0 ],
              "BB-SC ref": [ 0, 1 ],
              "BB-SC angle": 100,
              "vsites3": { 2: [ 0.35, 0.00, 1.053 ],
                           3: [ 0.00, 0.35, 1.053 ] },
              "exclusion": [ 2, 3 ],
              "int exclusion": False,
              "vsitesn": []
              },
          "martini3": {
              "beads": [ "BB", "SC1"],
              "mass": [ None, None ],
              "FF name": [ "P2", "SQ5n"],
              "coords": np.array([[ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00],
                                  [ 3.52791575e+00,  3.98701858e-08, -3.49091158e-08]]),
              "charge": [ 0, -1],
              "BB bonds": [],
              "SC bonds": np.array([[ 0, 1, 0.352, 7500 ]]),
              "constraints": [],
              "BB angles": [],
              "SC angles": [],
              "BB dihedral": [],
              "SC dihedral": [],
              "BB ref": [ 0 ],
              "BB-SC ref": [ 0, 1 ],
              "BB-SC angle": 100,
              "vsites3": {  },
              "exclusion": [ ],
              "int exclusion": False,
              "vsitesn": []
              },
          "marstini3": {
              "beads": [ "BB", "SC1", "H", "O" ],
              "mass": [ None, None, 0.00, 0.00 ],
              "FF name": [ "CPBB", "SQ5n", "CPNH", "CPCO" ],
              "coords": np.array([[ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00],
                                  [ 3.52791575e+00,  3.98701858e-08, -3.49091158e-08],
                                  [ 0.00000000e+00,  3.50000000e-01,  1.05300000e+00],
                                  [ 0.00000000e+00, -3.50000000e-01,  1.05300000e+00]]),
              "charge": [ 0, -1, 0.4, -0.4 ],
              "BB bonds": [],
              "SC bonds": np.array([[ 0, 1, 0.352, 7500 ]]),
              "constraints": [],
              "BB angles": [],
              "SC angles": [],
              "BB dihedral": [],
              "SC dihedral": [],
              "BB ref": [ 0 ],
              "BB-SC ref": [ 0, 1 ],
              "BB-SC angle": 100,
              "vsites3": { 2: [ 0.35, 0.00, 1.053 ],
                           3: [ 0.00, 0.35, 1.053 ] },
              "exclusion": [ 2, 3 ],
              "int exclusion": False,
              "vsitesn": []
              },
          "martini22p": {
              
              } 
        },
        
    # Glutamic acid
    "E" :
        { "name": "GLU",
          "AA" : {
              "beads": [ "N", "H", "CA", "CB", "CG", "CD", "OE1", "OE2", "C", "O" ],
              "coords": np.array([[ -7.070e-01, -1.225e+00,  5.000e-01 ],
                                  [ -7.070e-01, -1.225e+00,  1.500e+00 ],
                                  [  0.000e+00,  0.000e+00,  0.000e+00 ],
                                  [  1.520e+00,  0.000e+00,  5.070e-01 ],
                                  [  2.135e+00,  1.390e+00,  5.070e-01 ],
                                  [  3.654e+00,  1.376e+00,  5.020e-01 ],
                                  [  4.251e+00,  2.830e-01,  6.180e-01 ],
                                  [  4.252e+00,  2.467e+00,  3.820e-01 ],
                                  [ -7.160e-01,  1.241e+00,  5.070e-01 ],
                                  [ -7.160e-01,  1.241e+00,  1.737e+00 ]]),
              "BB ref": [ 0, 2, 8, 9 ],
              "SC ref": [ 3 ],
              },
          "martini22":  {
              "beads": [ "BB", "SC1" ],
              "mass": [ None, None ],
              "FF name": [ "P5", "Qa" ],
              "coords": np.array([[ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00],
                                  [ 4.35432938e+00,  1.25751584e-07,  2.85747300e-07]]),
              "charge": [ 0, -1 ],
              "BB bonds": [],
              "SC bonds": np.array([[ 0, 1, 0.40, 5000 ]]),
              "constraints": [],
              "BB angles": [],
              "SC angles": [],
              "BB dihedral": [],
              "SC dihedral": [],
              "BB ref": [ 0 ],
              "BB-SC ref": [ 0, 1 ],
              "BB-SC angle": 100,
              "vsites3": { },
              "exclusion": [ ],
              "int exclusion": False,
              "vsitesn": []
              },
          "marstini22": {
              "beads": [ "BB", "SC1", "H", "O" ],
              "mass": [ None, None, 0.00, 0.00 ],
              "FF name": [ "CPBB", "Qa", "CPNH", "CPCO" ],
              "coords": np.array([[ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00],
                                  [ 4.35432938e+00,  1.25751584e-07,  2.85747300e-07],
                                  [ 0.00000000e+00,  3.50000000e-01,  1.05300000e+00],
                                  [ 0.00000000e+00, -3.50000000e-01,  1.05300000e+00]]),
              "charge": [ 0, -1, 0.4, -0.4 ],
              "BB bonds": [],
              "SC bonds": np.array([[ 0, 1, 0.40, 5000 ]]),
              "constraints": [],
              "BB angles": [],
              "SC angles": [],
              "BB dihedral": [],
              "SC dihedral": [],
              "BB ref": [ 0 ],
              "BB-SC ref": [ 0, 1 ],
              "BB-SC angle": 100,
              "vsites3": { 2: [ 0.35, 0.00, 1.053 ],
                           3: [ 0.00, 0.35, 1.053 ] },
              "exclusion": [ 2, 3 ],
              "int exclusion": False,
              "vsitesn": []
              },
          "martini3": {
              "beads": [ "BB", "SC1"],
              "mass": [ None, None ],
              "FF name": [ "P2", "Q5n"],
              "coords": np.array([[ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00],
                                  [ 4.35432938e+00,  1.25751584e-07,  2.85747300e-07]]),
              "charge": [ 0, -1 ],
              "BB bonds": [],
              "SC bonds": np.array([[ 0, 1, 0.40, 5000 ]]),
              "constraints": [],
              "BB angles": [],
              "SC angles": [],
              "BB dihedral": [],
              "SC dihedral": [],
              "BB ref": [ 0 ],
              "BB-SC ref": [ 0, 1 ],
              "BB-SC angle": 100,
              "vsites3": { },
              "exclusion": [ ],
              "int exclusion": False,
              "vsitesn": []
              },
          "marstini3": {
              "beads": [ "BB", "SC1", "H", "O" ],
              "mass": [ None, None, 0.00, 0.00 ],
              "FF name": [ "CPBB", "Q5n", "CPNH", "CPCO" ],
              "coords": np.array([[ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00],
                                  [ 4.35432938e+00,  1.25751584e-07,  2.85747300e-07],
                                  [ 0.00000000e+00,  3.50000000e-01,  1.05300000e+00],
                                  [ 0.00000000e+00, -3.50000000e-01,  1.05300000e+00]]),
              "charge": [ 0, -1, 0.4, -0.4 ],
              "BB bonds": [],
              "SC bonds": np.array([[ 0, 1, 0.40, 5000 ]]),
              "constraints": [],
              "BB angles": [],
              "SC angles": [],
              "BB dihedral": [],
              "SC dihedral": [],
              "BB ref": [ 0 ],
              "BB-SC ref": [ 0, 1 ],
              "BB-SC angle": 100,
              "vsites3": { 2: [ 0.35, 0.00, 1.053 ],
                           3: [ 0.00, 0.35, 1.053 ] },
              "exclusion": [ 2, 3 ],
              "int exclusion": False,
              "vsitesn": []
              },
          "martini22p": {
              
              } 
        },
        
    # Serine
    "S" :
        { "name": "SER",
          "AA" : {
              "beads": [ "N", "H", "CA", "CB", "OG", "HG", "C", "O" ],
              "coords": np.array([[ -7.070e-01, -1.225e+00,  5.000e-01 ],
                                  [ -7.070e-01, -1.225e+00,  1.500e+00 ],
                                  [  0.000e+00,  0.000e+00,  0.000e+00 ],
                                  [  1.520e+00,  0.000e+00,  5.070e-01 ],
                                  [  2.023e+00,  1.326e+00,  5.070e-01 ],
                                  [  3.022e+00,  1.303e+00,  5.070e-01 ],
                                  [ -7.160e-01,  1.241e+00,  5.070e-01 ],
                                  [ -7.160e-01,  1.241e+00,  1.737e+00 ]]),
              "BB ref": [ 0, 2, 6, 7 ],
              "SC ref": [ 3 ],
              },
          "martini22":  {
              "beads": [ "BB", "SC1" ],
              "mass": [ None, None ],
              "FF name": [ "P5", "P1" ],
              "coords": np.array([[ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00],
                                  [ 2.88263697e+00,  5.03615305e-08,  7.97614532e-08]]),
              "charge": [ 0, 0 ],
              "BB bonds": [],
              "SC bonds": np.array([[ 0, 1, 0.25, 7500 ]]),
              "constraints": [],
              "BB angles": [],
              "SC angles": [],
              "BB dihedral": [],
              "SC dihedral": [],
              "BB ref": [ 0 ],
              "BB-SC ref": [ 0, 1 ],
              "BB-SC angle": 100,
              "vsites3": { },
              "exclusion": [ ],
              "int exclusion": False,
              "vsitesn": []
              },
          "marstini22": {
              "beads": [ "BB", "SC1", "H", "O" ],
              "mass": [ None, None, 0.00, 0.00 ],
              "FF name": [ "CPBB", "P1", "CPNH", "CPCO" ],
              "coords": np.array([[ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00],
                                  [ 2.88263697e+00,  5.03615305e-08,  7.97614532e-08],
                                  [ 0.00000000e+00,  3.50000000e-01,  1.05300000e+00],
                                  [ 0.00000000e+00, -3.50000000e-01,  1.05300000e+00]]),
              "charge": [ 0, 0, 0.4, -0.4 ],
              "BB bonds": [],
              "SC bonds": np.array([[ 0, 1, 0.25, 7500 ]]),
              "constraints": [],
              "BB angles": [],
              "SC angles": [],
              "BB dihedral": [],
              "SC dihedral": [],
              "BB ref": [ 0 ],
              "BB-SC ref": [ 0, 1 ],
              "BB-SC angle": 100,
              "vsites3": { 2: [ 0.35, 0.00, 1.053 ],
                           3: [ 0.00, 0.35, 1.053 ] },
              "exclusion": [ 2, 3 ],
              "int exclusion": False,
              "vsitesn": []
              },
          "martini3": {
              "beads": [ "BB", "SC1"],
              "mass": [ None, None ],
              "FF name": [ "P2", "TP1" ],
              "coords": np.array([[ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00],
                                  [ 2.88263697e+00,  5.03615305e-08,  7.97614532e-08]]),
              "charge": [ 0, 0 ],
              "BB bonds": [],
              "SC bonds": np.array([[ 0, 1, 0.287, 7500 ]]),
              "constraints": [],
              "BB angles": [],
              "SC angles": [],
              "BB dihedral": [],
              "SC dihedral": [],
              "BB ref": [ 0 ],
              "BB-SC ref": [ 0, 1 ],
              "BB-SC angle": 100,
              "vsites3": { },
              "exclusion": [  ],
              "int exclusion": False,
              "vsitesn": []
              },
          "marstini3": {
              "beads": [ "BB", "SC1", "H", "O" ],
              "mass": [ None, None, 0.00, 0.00 ],
              "FF name": [ "CPBB", "TP1", "CPNH", "CPCO" ],
              "coords": np.array([[ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00],
                                  [ 2.88263697e+00,  5.03615305e-08,  7.97614532e-08],
                                  [ 0.00000000e+00,  3.50000000e-01,  1.05300000e+00],
                                  [ 0.00000000e+00, -3.50000000e-01,  1.05300000e+00]]),
              "charge": [ 0, 0, 0.4, -0.4 ],
              "BB bonds": [],
              "SC bonds": np.array([[ 0, 1, 0.287, 7500 ]]),
              "constraints": [],
              "BB angles": [],
              "SC angles": [],
              "BB dihedral": [],
              "SC dihedral": [],
              "BB ref": [ 0 ],
              "BB-SC ref": [ 0, 1 ],
              "BB-SC angle": 100,
              "vsites3": { 2: [ 0.35, 0.00, 1.053 ],
                           3: [ 0.00, 0.35, 1.053 ] },
              "exclusion": [ 2, 3 ],
              "int exclusion": False,
              "vsitesn": []
              },
          "martini22p": {
              
              } 
        },
        
    # Threonine
    "T" :
        { "name": "THR",
          "AA" : {
              "beads": [ "N", "H", "CA", "CB", "OG1", "HG1", "CG2", "C", "O" ],
              "coords": np.array([[ -7.071e-01, -1.225e+00,  5.000e-01 ],
                                  [ -7.071e-01, -1.225e+00,  1.500e+00 ],
                                  [  0.000e+00,  0.000e+00,  0.000e+00 ],
                                  [  1.520e+00,  2.672e-09,  5.067e-01 ],
                                  [  1.989e+00,  1.271e+00,  5.083e-02 ],
                                  [  2.990e+00,  1.276e+00,  4.888e-02 ],
                                  [  2.070e+00, -1.094e+00, -4.108e-01 ],
                                  [ -7.165e-01,  1.241e+00,  5.067e-01 ],
                                  [ -7.165e-01,  1.241e+00,  1.737e+00 ]]),
              "BB ref": [ 0, 2, 7, 8 ],
              "SC ref": [ 3 ],
              },
          "martini22":  {
              "beads": [ "BB", "SC1" ],
              "mass": [ None, None ],
              "FF name": [ "P5", "P1" ],
              "coords": np.array([[ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00],
                                  [ 2.64803708e+00,  4.72833052e-08, -1.27962501e-07]]),
              "charge": [ 0, 0 ],
              "BB bonds": [],
              "SC bonds": [],
              "constraints": np.array([[ 0, 1, 0.26 ]]),
              "BB angles": [],
              "SC angles": [],
              "BB dihedral": [],
              "SC dihedral": [],
              "BB ref": [ 0 ],
              "BB-SC ref": [ 0, 1 ],
              "BB-SC angle": 100,
              "vsites3": { },
              "exclusion": [ ],
              "int exclusion": False,
              "vsitesn": []
              },
          "marstini22": {
              "beads": [ "BB", "SC1", "H", "O" ],
              "mass": [ None, None, 0.00, 0.00 ],
              "FF name": [ "CPBB", "P1", "CPNH", "CPCO" ],
              "coords": np.array([[ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00],
                                  [ 2.64803708e+00,  4.72833052e-08, -1.27962501e-07],
                                  [ 0.00000000e+00,  3.50000000e-01,  1.05300000e+00],
                                  [ 0.00000000e+00, -3.50000000e-01,  1.05300000e+00]]),
              "charge": [ 0, 0, 0.4, -0.4 ],
              "BB bonds": [],
              "SC bonds": [],
              "constraints": np.array([[ 0, 1, 0.26 ]]),
              "BB angles": [],
              "SC angles": [],
              "BB dihedral": [],
              "SC dihedral": [],
              "BB ref": [ 0 ],
              "BB-SC ref": [ 0, 1 ],
              "BB-SC angle": 100,
              "vsites3": { 2: [ 0.35, 0.00, 1.053 ],
                           3: [ 0.00, 0.35, 1.053 ] },
              "exclusion": [ 2, 3 ],
              "int exclusion": False,
              "vsitesn": []
              },
          "martini3": {
              "beads": [ "BB", "SC1" ],
              "mass": [ None, None ],
              "FF name": [ "P2", "SP1"],
              "coords": np.array([[ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00],
                                  [ 2.64803708e+00,  4.72833052e-08, -1.27962501e-07]]),
              "charge": [ 0, 0 ],
              "BB bonds": [],
              "SC bonds": [],
              "constraints": np.array([[ 0, 1, 0.305 ]]),
              "BB angles": [],
              "SC angles": [],
              "BB dihedral": [],
              "SC dihedral": [],
              "BB ref": [ 0 ],
              "BB-SC ref": [ 0, 1 ],
              "BB-SC angle": 100,
              "vsites3": { },
              "exclusion": [  ],
              "int exclusion": False,
              "vsitesn": []
              },
          "marstini3": {
              "beads": [ "BB", "SC1", "H", "O" ],
              "mass": [ None, None, 0.00, 0.00 ],
              "FF name": [ "CPBB", "SP1", "CPNH", "CPCO" ],
              "coords": np.array([[ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00],
                                  [ 2.64803708e+00,  4.72833052e-08, -1.27962501e-07],
                                  [ 0.00000000e+00,  3.50000000e-01,  1.05300000e+00],
                                  [ 0.00000000e+00, -3.50000000e-01,  1.05300000e+00]]),
              "charge": [ 0, 0, 0.4, -0.4 ],
              "BB bonds": [],
              "SC bonds": [],
              "constraints": np.array([[ 0, 1, 0.305 ]]),
              "BB angles": [],
              "SC angles": [],
              "BB dihedral": [],
              "SC dihedral": [],
              "BB ref": [ 0 ],
              "BB-SC ref": [ 0, 1 ],
              "BB-SC angle": 100,
              "vsites3": { 2: [ 0.35, 0.00, 1.053 ],
                           3: [ 0.00, 0.35, 1.053 ] },
              "exclusion": [ 2, 3 ],
              "int exclusion": False,
              "vsitesn": []
              },
          "martini22p": {

              } 
        },
        
    # Asparagine
    "N" :
        { "name": "ASN",
          "AA" : {
              "beads": [ "N", "H", "CA", "CB", "CG", "OD1", "ND2", "HD21", "HD22", "C", "O" ],
              "coords": np.array([[ -7.070e-01, -1.225e+00,  5.000e-01 ],
                                  [ -7.070e-01, -1.225e+00,  1.500e+00 ],
                                  [  0.000e+00,  0.000e+00,  0.000e+00 ],
                                  [  1.520e+00,  0.000e+00,  5.070e-01 ],
                                  [  2.105e+00,  1.403e+00,  5.070e-01 ],
                                  [  1.373e+00,  2.390e+00,  5.070e-01 ],
                                  [  3.432e+00,  1.491e+00,  5.060e-01 ],
                                  [  3.989e+00,  6.600e-01,  5.050e-01 ],
                                  [  3.874e+00,  2.389e+00,  5.060e-01 ],
                                  [ -7.160e-01,  1.241e+00,  5.070e-01 ],
                                  [ -7.160e-01,  1.241e+00,  1.737e+00 ]]),
              "BB ref": [ 0, 2, 9, 10 ],
              "SC ref": [ 3 ],
              },
          "martini22":  {
              "beads": [ "BB", "SC1" ],
              "mass": [ None, None ],
              "FF name": [ "P5", "P5" ],
              "coords": np.array([[ 0.00000000e+00, 0.00000000e+00, 0.00000000e+00],
                                  [ 3.49509616e+00, 6.30210431e-08, 5.93628864e-08]]),
              "charge": [ 0, 0 ],
              "BB bonds": [],
              "SC bonds": np.array([[ 0, 1, 0.32, 5000 ]]),
              "constraints": [],
              "BB angles": [],
              "SC angles": [],
              "BB dihedral": [],
              "SC dihedral": [],
              "BB ref": [ 0 ],
              "BB-SC ref": [ 0, 1 ],
              "BB-SC angle": 100,
              "vsites3": { },
              "exclusion": [ ],
              "int exclusion": False,
              "vsitesn": []
              },
          "marstini22": {
              "beads": [ "BB", "SC1", "H", "O" ],
              "mass": [ None, None, 0.00, 0.00 ],
              "FF name": [ "CPBB", "P5", "CPNH", "CPCO" ],
              "coords": np.array([[ 0.00000000e+00, 0.00000000e+00, 0.00000000e+00],
                                  [ 3.49509616e+00, 6.30210431e-08, 5.93628864e-08],
                                  [ 0.00000000e+00,  3.50000000e-01,  1.05300000e+00],
                                  [ 0.00000000e+00, -3.50000000e-01,  1.05300000e+00]]),
              "charge": [ 0, 0, 0.4, -0.4 ],
              "BB bonds": [],
              "SC bonds": np.array([[ 0, 1, 0.32, 5000 ]]),
              "constraints": [],
              "BB angles": [],
              "SC angles": [],
              "BB dihedral": [],
              "SC dihedral": [],
              "BB ref": [ 0 ],
              "BB-SC ref": [ 0, 1 ],
              "BB-SC angle": 100,
              "vsites3": { 2: [ 0.35, 0.00, 1.053 ],
                           3: [ 0.00, 0.35, 1.053 ] },
              "exclusion": [ 2, 3 ],
              "int exclusion": False,
              "vsitesn": []
              },
          "martini3": {
              "beads": [ "BB", "SC1" ],
              "mass": [ None, None ],
              "FF name": [ "P2", "SP5" ],
              "coords": np.array([[ 0.00000000e+00, 0.00000000e+00, 0.00000000e+00],
                                  [ 3.49509616e+00, 6.30210431e-08, 5.93628864e-08]]),
              "charge": [ 0, 0 ],
              "BB bonds": [],
              "SC bonds": np.array([[ 0, 1, 0.352, 5000 ]]),
              "constraints": [],
              "BB angles": [],
              "SC angles": [],
              "BB dihedral": [],
              "SC dihedral": [],
              "BB ref": [ 0 ],
              "BB-SC ref": [ 0, 1 ],
              "BB-SC angle": 100,
              "vsites3": { },
              "exclusion": [ ],
              "int exclusion": False,
              "vsitesn": []
              },
          "marstini3": {
              "beads": [ "BB", "SC1", "H", "O" ],
              "mass": [ None, None, 0.00, 0.00 ],
              "FF name": [ "CPBB", "SP5", "CPNH", "CPCO" ],
              "coords": np.array([[ 0.00000000e+00, 0.00000000e+00, 0.00000000e+00],
                                  [ 3.49509616e+00, 6.30210431e-08, 5.93628864e-08],
                                  [ 0.00000000e+00,  3.50000000e-01,  1.05300000e+00],
                                  [ 0.00000000e+00, -3.50000000e-01,  1.05300000e+00]]),
              "charge": [ 0, 0, 0.4, -0.4 ],
              "BB bonds": [],
              "SC bonds": np.array([[ 0, 1, 0.352, 5000 ]]),
              "constraints": [],
              "BB angles": [],
              "SC angles": [],
              "BB dihedral": [],
              "SC dihedral": [],
              "BB ref": [ 0 ],
              "BB-SC ref": [ 0, 1 ],
              "BB-SC angle": 100,
              "vsites3": { 2: [ 0.35, 0.00, 1.053 ],
                           3: [ 0.00, 0.35, 1.053 ] },
              "exclusion": [ 2, 3 ],
              "int exclusion": False,
              "vsitesn": []
              },
          "martini22p": {
              
              } 
        },
        
    # Glutamine
    "Q" :
        { "name": "GLN",
          "AA" : {
              "beads": [ "N", "H", "CA", "CB", "CG", "CD", "OE1", "NE2", "HE21", "HE22", "C", "O" ],
              "coords": np.array([[ -7.070e-01, -1.225e+00,  5.000e-01 ],
                                  [ -7.070e-01, -1.225e+00,  1.500e+00 ],
                                  [  0.000e+00,  0.000e+00,  0.000e+00 ],
                                  [  1.519e+00, -0.000e+00,  5.070e-01 ],
                                  [  2.132e+00,  1.391e+00,  5.070e-01 ],
                                  [  3.652e+00,  1.365e+00,  5.070e-01 ],
                                  [  4.306e+00,  2.417e+00,  5.080e-01 ],
                                  [  4.225e+00,  1.640e-01,  5.070e-01 ],
                                  [  3.658e+00, -6.600e-01,  5.070e-01 ],
                                  [  5.221e+00,  8.600e-02,  5.080e-01 ],
                                  [ -7.160e-01,  1.241e+00,  5.070e-01 ],
                                  [ -7.160e-01,  1.241e+00,  1.737e+00 ]]),
              "BB ref": [ 0, 2, 10, 11 ],
              "SC ref": [ 3 ],
              },
          "martini22":  {
              "beads": [ "BB", "SC1" ],
              "mass": [ None, None ],
              "FF name": [ "P5", "Qa" ],
              "coords": np.array([[ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00],
                                  [ 3.82478721e+00, -5.06190360e-07,  1.65297370e-07]]),
              "charge": [ 0, 0 ],
              "BB bonds": [],
              "SC bonds": np.array([[ 0, 1, 0.40, 5000 ]]),
              "constraints": [],
              "BB angles": [],
              "SC angles": [],
              "BB dihedral": [],
              "SC dihedral": [],
              "BB ref": [ 0 ],
              "BB-SC ref": [ 0, 1 ],
              "BB-SC angle": 100,
              "vsites3": { },
              "exclusion": [ ],
              "int exclusion": False,
              "vsitesn": []
              },
          "marstini22": {
              "beads": [ "BB", "SC1", "H", "O" ],
              "mass": [ None, None, 0.00, 0.00 ],
              "FF name": [ "CPBB", "Qa", "CPNH", "CPCO" ],
              "coords": np.array([[ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00],
                                  [ 3.82478721e+00, -5.06190360e-07,  1.65297370e-07],
                                  [ 0.00000000e+00,  3.50000000e-01,  1.05300000e+00],
                                  [ 0.00000000e+00, -3.50000000e-01,  1.05300000e+00]]),
              "charge": [ 0, 0, 0.4, -0.4 ],
              "BB bonds": [],
              "SC bonds": np.array([[ 0, 1, 0.40, 5000 ]]),
              "constraints": [],
              "BB angles": [],
              "SC angles": [],
              "BB dihedral": [],
              "SC dihedral": [],
              "BB ref": [ 0 ],
              "BB-SC ref": [ 0, 1 ],
              "BB-SC angle": 100,
              "vsites3": { 2: [ 0.35, 0.00, 1.053 ],
                           3: [ 0.00, 0.35, 1.053 ] },
              "exclusion": [ 2, 3 ],
              "int exclusion": False,
              "vsitesn": []
              },
          "martini3": {
              "beads": [ "BB", "SC1"],
              "mass": [ None, None ],
              "FF name": [ "P2", "P5" ],
              "coords": np.array([[ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00],
                                  [ 3.82478721e+00, -5.06190360e-07,  1.65297370e-07]]),
              "charge": [ 0, 0 ],
              "BB bonds": [],
              "SC bonds": np.array([[ 0, 1, 0.40, 5000 ]]),
              "constraints": [],
              "BB angles": [],
              "SC angles": [],
              "BB dihedral": [],
              "SC dihedral": [],
              "BB ref": [ 0 ],
              "BB-SC ref": [ 0, 1 ],
              "BB-SC angle": 100,
              "vsites3": { },
              "exclusion": [ ],
              "int exclusion": False,
              "vsitesn": []
              },
          "marstini3": {
              "beads": [ "BB", "SC1", "H", "O" ],
              "mass": [ None, None, 0.00, 0.00 ],
              "FF name": [ "CPBB", "P5", "CPNH", "CPCO" ],
              "coords": np.array([[ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00],
                                  [ 3.82478721e+00, -5.06190360e-07,  1.65297370e-07],
                                  [ 0.00000000e+00,  3.50000000e-01,  1.05300000e+00],
                                  [ 0.00000000e+00, -3.50000000e-01,  1.05300000e+00]]),
              "charge": [ 0, 0, 0.4, -0.4 ],
              "BB bonds": [],
              "SC bonds": np.array([[ 0, 1, 0.40, 5000 ]]),
              "constraints": [],
              "BB angles": [],
              "SC angles": [],
              "BB dihedral": [],
              "SC dihedral": [],
              "BB ref": [ 0 ],
              "BB-SC ref": [ 0, 1 ],
              "BB-SC angle": 100,
              "vsites3": { 2: [ 0.35, 0.00, 1.053 ],
                           3: [ 0.00, 0.35, 1.053 ] },
              "exclusion": [ 2, 3 ],
              "int exclusion": False,
              "vsitesn": []
              },
          "martini22p": {
              
              } 
        },
        
    # Cysteine
    "C" :
        { "name": "CYS",
          "AA" : {
              "beads": [ "N", "H", "CA", "CB", "SG", "HG", "C", "O" ],
              "coords": np.array([[ -7.070e-01, -1.225e+00,  5.000e-01 ],
                                  [ -7.070e-01, -1.225e+00,  1.500e+00 ],
                                  [  0.000e+00,  0.000e+00,  0.000e+00 ],
                                  [  1.520e+00,  0.000e+00,  5.070e-01 ],
                                  [  2.250e+00,  1.654e+00,  5.070e-01 ],
                                  [  3.247e+00,  1.579e+00,  5.070e-01 ],
                                  [ -7.160e-01,  1.241e+00,  5.070e-01 ],
                                  [ -7.160e-01,  1.241e+00,  1.737e+00 ]]),
              "BB ref": [ 0, 2, 6, 7 ],
              "SC ref": [ 3 ],
              },
          "martini22":  {
              "beads": [ "BB", "SC1" ],
              "mass": [ None, None ],
              "FF name": [ "P5", "C5" ],
              "coords": np.array([[ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00],
                                  [ 3.32928011e+00,  2.72712204e-08,  2.98065752e-08]]),
              "charge": [ 0, 0 ],
              "BB bonds": [],
              "SC bonds": np.array([[ 0, 1, 0.31, 7500 ]]),
              "constraints": [],
              "BB angles": [],
              "SC angles": [],
              "BB dihedral": [],
              "SC dihedral": [],
              "BB ref": [ 0 ],
              "BB-SC ref": [ 0, 1 ],
              "BB-SC angle": 100,
              "vsites3": { },
              "exclusion": [ ],
              "int exclusion": False,
              "vsitesn": []
              },
          "marstini22": {
              "beads": [ "BB", "SC1" , "H", "O" ],
              "mass": [ None, None, 0.00, 0.00 ],
              "FF name": [ "CPBB", "C5", "CPNH", "CPCO" ],
              "coords": np.array([[ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00],
                                  [ 3.32928011e+00,  2.72712204e-08,  2.98065752e-08],
                                  [ 0.00000000e+00,  3.50000000e-01,  1.05300000e+00],
                                  [ 0.00000000e+00, -3.50000000e-01,  1.05300000e+00]]),
              "charge": [ 0, 0, 0.4, -0.4 ],
              "BB bonds": [],
              "SC bonds": np.array([[ 0, 1, 0.31, 7500 ]]),
              "constraints": [],
              "BB angles": [],
              "SC angles": [],
              "BB dihedral": [],
              "SC dihedral": [],
              "BB ref": [ 0 ],
              "BB-SC ref": [ 0, 1 ],
              "BB-SC angle": 100,
              "vsites3": { 2: [ 0.35, 0.00, 1.053 ],
                           3: [ 0.00, 0.35, 1.053 ] },
              "exclusion": [ 2, 3 ],
              "int exclusion": False,
              "vsitesn": []
              },
          "martini3": {
              "beads": [ "BB", "SC1"],
              "mass": [ None, None ],
              "FF name": [ "P2", "TC6" ],
              "coords": np.array([[ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00],
                                  [ 3.32928011e+00,  2.72712204e-08,  2.98065752e-08]]),
              "charge": [ 0, 0 ],
              "BB bonds": [],
              "SC bonds": np.array([[ 0, 1, 0.341, 7500 ]]),
              "constraints": [],
              "BB angles": [],
              "SC angles": [],
              "BB dihedral": [],
              "SC dihedral": [],
              "BB ref": [ 0 ],
              "BB-SC ref": [ 0, 1 ],
              "BB-SC angle": 100,
              "vsites3": { },
              "exclusion": [ ],
              "int exclusion": False,
              "vsitesn": []
              },
          "marstini3": {
              "beads": [ "BB", "SC1" , "H", "O" ],
              "mass": [ None, None, 0.00, 0.00 ],
              "FF name": [ "CPBB", "TC6", "CPNH", "CPCO" ],
              "coords": np.array([[ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00],
                                  [ 3.32928011e+00,  2.72712204e-08,  2.98065752e-08],
                                  [ 0.00000000e+00,  3.50000000e-01,  1.05300000e+00],
                                  [ 0.00000000e+00, -3.50000000e-01,  1.05300000e+00]]),
              "charge": [ 0, 0, 0.4, -0.4 ],
              "BB bonds": [],
              "SC bonds": np.array([[ 0, 1, 0.341, 7500 ]]),
              "constraints": [],
              "BB angles": [],
              "SC angles": [],
              "BB dihedral": [],
              "SC dihedral": [],
              "BB ref": [ 0 ],
              "BB-SC ref": [ 0, 1 ],
              "BB-SC angle": 100,
              "vsites3": { 2: [ 0.35, 0.00, 1.053 ],
                           3: [ 0.00, 0.35, 1.053 ] },
              "exclusion": [ 2, 3 ],
              "int exclusion": False,
              "vsitesn": []
              },
          "martini22p": {
              
              }
        },
        
    # Glycine
    "G" :
        { "name": "GLY",
          "AA": {
              "beads": [ "N", "H", "CA", "C", "O" ],
              "coords": np.array([[ -7.071e-01, -1.225e+00,  5.000e-01 ],
                                  [ -7.071e-01, -1.225e+00,  1.500e+00 ],
                                  [  0.000e+00,  0.000e+00,  0.000e+00 ],
                                  [ -7.165e-01,  1.241e+00,  5.067e-01 ],
                                  [ -7.165e-01,  1.241e+00,  1.737e+00 ]]),
              "BB ref": [ 0, 2, 5, 6 ],
              "SC ref": [  ],
              },
          "martini22":  {
              "beads": [ "BB" ],
              "mass": [ None ],
              "FF name": [ "P5" ],
              "coords": np.array([[ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00]]),
              "charge": [ 0 ],
              "BB bonds": [],
              "SC bonds": [],
              "constraints": [],
              "BB angles": [],
              "SC angles": [],
              "BB dihedral": [],
              "SC dihedral": [],
              "BB ref": [ 0 ],
              "BB-SC ref": [],
              "BB-SC angle": 100,
              "vsites3": { },
              "exclusion": [ ],
              "int exclusion": False,
              "vsitesn": []
              },
          "marstini22": {
              "beads": [ "BB", "H", "O" ],
              "mass": [ None, 0.00, 0.00 ],
              "FF name": [ "CPBB", "CPNH", "CPCO" ],
              "coords": np.array([[ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00],
                                  [ 0.00000000e+00,  3.50000000e-01,  1.05300000e+00],
                                  [ 0.00000000e+00, -3.50000000e-01,  1.05300000e+00]]),
              "charge": [ 0, 0.4, -0.4 ],
              "BB bonds": [],
              "SC bonds": [],
              "constraints": [],
              "BB angles": [],
              "SC angles": [],
              "BB dihedral": [],
              "SC dihedral": [],
              "BB ref": [ 0 ],
              "BB-SC ref": [],
              "BB-SC angle": 100,
              "vsites3": { 1: [ 0.35, 0.00, 1.053 ],
                           2: [ 0.00, 0.35, 1.053 ] },
              "exclusion": [ 1, 2 ],
              "int exclusion": False,
              "vsitesn": []
              },
          "martini3": {
              "beads": [ "BB" ],
              "mass": [ None ],
              "FF name": [ "SP1" ],
              "coords": np.array([[ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00]]),
              "charge": [ 0 ],
              "BB bonds": [],
              "SC bonds": [],
              "constraints": [],
              "BB angles": [],
              "SC angles": [],
              "BB dihedral": [],
              "SC dihedral": [],
              "BB ref": [ 0 ],
              "BB-SC ref": [],
              "BB-SC angle": 100,
              "vsites3": {},
              "exclusion": [ ],
              "int exclusion": False,
              "vsitesn": []
              },
          "marstini3": {
              "beads": [ "BB", "H", "O" ],
              "mass": [ None, 0.00, 0.00 ],
              "FF name": [ "CPBB", "CPNH", "CPCO" ],
              "coords": np.array([[ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00],
                                  [ 0.00000000e+00,  3.50000000e-01,  1.05300000e+00],
                                  [ 0.00000000e+00, -3.50000000e-01,  1.05300000e+00]]),
              "charge": [ 0, 0.4, -0.4 ],
              "BB bonds": [],
              "SC bonds": [],
              "constraints": [],
              "BB angles": [],
              "SC angles": [],
              "BB dihedral": [],
              "SC dihedral": [],
              "BB ref": [ 0 ],
              "BB-SC ref": [],
              "BB-SC angle": 100,
              "vsites3": { 1: [ 0.35, 0.00, 1.053 ],
                           2: [ 0.00, 0.35, 1.053 ] },
              "exclusion": [ 1, 2 ],
              "int exclusion": False,
              "vsitesn": []
              },
          "martini22p": {
              
              } 
        },
        
    # Proline
    "P" :
        { "name": "PRO",
          "AA" : {
              "beads": [ "N", "CA", "CB", "CG", "CD", "C", "O" ],
              "coords": np.array([[ -7.690e-01, -1.200e+00,  3.170e-01 ],
                                  [ -1.500e-01, -1.929e+00,  1.446e+00 ],
                                  [  0.000e+00,  0.000e+00,  0.000e+00 ],
                                  [  1.309e+00, -1.920e-01,  7.480e-01 ],
                                  [  9.340e-01, -1.014e+00,  1.933e+00 ],
                                  [ -6.750e-01,  1.275e+00,  4.790e-01 ],
                                  [ -8.320e-01,  1.485e+00,  1.680e+00 ]]),
              "BB ref": [ 0, 1, 5, 6 ],
              "SC ref": [ 2 ],
              },
          "martini22":  {
              "beads": [ "BB", "SC1" ],
              "mass": [ None, None ],
              "FF name": [ "P4", "C3" ],
              "coords": np.array([[ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00],
                                  [ 2.26267503e+00, -5.29634170e-09,  7.06946801e-09]]),
              "charge": [ 0, 0 ],
              "BB bonds": [],
              "SC bonds": np.array([[ 0, 1, 0.3, 7500 ]]),
              "constraints": [],
              "BB angles": [],
              "SC angles": [],
              "BB dihedral": [],
              "SC dihedral": [],
              "BB ref": [ 0 ],
              "BB-SC ref": [],
              "BB-SC angle": 100,
              "vsites3": { },
              "exclusion": [ ],
              "int exclusion": False,
              "vsitesn": []
              },
          "marstini22": {
              "beads": [ "BB", "SC1" , "H", "O" ],
              "mass": [ None, None, 0.00, 0.00 ],
              "FF name": [ "CPBB", "C3", "CPNH", "CPCO" ],
              "coords": np.array([[ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00],
                                  [ 2.26267503e+00, -5.29634170e-09,  7.06946801e-09],
                                  [ 0.00000000e+00,  3.50000000e-01,  1.05300000e+00],
                                  [ 0.00000000e+00, -3.50000000e-01,  1.05300000e+00]]),
              "charge": [ 0, 0, 0.4, -0.4 ],
              "BB bonds": [],
              "SC bonds": np.array([[ 0, 1, 0.3, 7500 ]]),
              "constraints": [],
              "BB angles": [],
              "SC angles": [],
              "BB dihedral": [],
              "SC dihedral": [],
              "BB ref": [ 0 ],
              "BB-SC ref": [],
              "BB-SC angle": 100,
              "vsites3": { 2: [ 0.35, 0.00, 1.053 ],
                           3: [ 0.00, 0.35, 1.053 ] },
              "exclusion": [ 2, 3 ],
              "int exclusion": False,
              "vsitesn": []
              },
          "martini3": {
              "beads": [ "BB", "SC1"],
              "mass": [ None, None ],
              "FF name": [ "P2", "C3" ],
              "coords": np.array([[ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00],
                                  [ 2.26267503e+00, -5.29634170e-09,  7.06946801e-09]]),
              "charge": [ 0, 0 ],
              "BB bonds": [],
              "SC bonds": np.array([[ 0, 1, 0.3, 7500 ]]),
              "constraints": [],
              "BB angles": [],
              "SC angles": [],
              "BB dihedral": [],
              "SC dihedral": [],
              "BB ref": [ 0 ],
              "BB-SC ref": [],
              "BB-SC angle": 100,
              "vsites3": { },
              "exclusion": [ ],
              "int exclusion": False,
              "vsitesn": []
              },
          "marstini3": {
              "beads": [ "BB", "SC1" , "H", "O" ],
              "mass": [ None, None, 0.00, 0.00 ],
              "FF name": [ "CPBB", "C3", "CPNH", "CPCO" ],
              "coords": np.array([[ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00],
                                  [ 2.26267503e+00, -5.29634170e-09,  7.06946801e-09],
                                  [ 0.00000000e+00,  3.50000000e-01,  1.05300000e+00],
                                  [ 0.00000000e+00, -3.50000000e-01,  1.05300000e+00]]),
              "charge": [ 0, 0, 0.4, -0.4 ],
              "BB bonds": [],
              "SC bonds": np.array([[ 0, 1, 0.3, 7500 ]]),
              "constraints": [],
              "BB angles": [],
              "SC angles": [],
              "BB dihedral": [],
              "SC dihedral": [],
              "BB ref": [ 0 ],
              "BB-SC ref": [],
              "BB-SC angle": 100,
              "vsites3": { 2: [ 0.35, 0.00, 1.053 ],
                           3: [ 0.00, 0.35, 1.053 ] },
              "exclusion": [ 2, 3 ],
              "int exclusion": False,
              "vsitesn": []
              },
          "martini22p": {
              
              } 
        },
        
    # Alanine
    "A" :
        { "name": "ALA",
          "AA" : {
              "beads": [ "N", "H", "CA", "CB", "C", "O" ],
              "coords": np.array([[ -7.071e-01, -1.225e+00,  5.000e-01 ],
                                  [ -7.071e-01, -1.225e+00,  1.500e+00 ],
                                  [  0.000e+00,  0.000e+00,  0.000e+00 ],
                                  [  1.520e+00,  1.653e-08,  5.067e-01 ],
                                  [ -7.165e-01,  1.241e+00,  5.067e-01 ],
                                  [ -7.165e-01,  1.241e+00,  1.737e+00 ]]),
              "BB ref": [ 0, 2, 4, 5 ],
              "SC ref": [ 3 ],
              },
          "martini22":  {
              "beads": [ "BB" ],
              "mass": [ None ],
              "FF name": [ "P4" ],
              "coords": np.array([[ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00]]),
              "charge": [ 0 ],
              "BB bonds": [],
              "SC bonds": [],
              "constraints": [],
              "BB angles": [],
              "SC angles": [],
              "BB dihedral": [],
              "SC dihedral": [],
              "BB ref": [ 0 ],
              "BB-SC ref": [],
              "BB-SC angle": 100,
              "vsites3": { },
              "exclusion": [ ],
              "int exclusion": False,
              "vsitesn": []
              },
          "marstini22": {
              "beads": [ "BB", "H", "O" ],
              "mass": [ None, 0.00, 0.00 ],
              "FF name": [ "CPBB", "CPNH", "CPCO" ],
              "coords": np.array([[ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00],
                                  [ 0.00000000e+00,  3.50000000e-01,  1.05300000e+00],
                                  [ 0.00000000e+00, -3.50000000e-01,  1.05300000e+00]]),
              "charge": [ 0, 0.4, -0.4 ],
              "BB bonds": [],
              "SC bonds": [],
              "constraints": [],
              "BB angles": [],
              "SC angles": [],
              "BB dihedral": [],
              "SC dihedral": [],
              "BB ref": [ 0 ],
              "BB-SC ref": [],
              "BB-SC angle": 100,
              "vsites3": { 1: [ 0.35, 0.00, 1.053 ],
                           2: [ 0.00, 0.35, 1.053 ] },
              "exclusion": [ 1, 2 ],
              "int exclusion": False,
              "vsitesn": []
              },
          "martini3": {
              "beads": [ "BB", "SC1"],
              "mass": [ None, None ],
              "FF name": [ "SP2", "TC3" ],
              "coords": np.array([[ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00],
                                  [ 2.83557758e+00,  5.34189583e-08,  3.07968850e-07]]),
              "charge": [ 0, 0 ],
              "BB bonds": [],
              "SC bonds": [],
              "constraints": np.array([[ 0, 1, 0.27]]),
              "BB angles": [],
              "SC angles": [],
              "BB dihedral": [],
              "SC dihedral": [],
              "BB ref": [ 0 ],
              "BB-SC ref": [ 0, 1 ],
              "BB-SC angle": 100,
              "vsites3": { },
              "exclusion": [],
              "int exclusion": False,
              "vsitesn": []
              },
          "marstini3": {
              "beads": [ "BB", "SC1", "H", "O" ],
              "mass": [ None, None, 0.00, 0.00 ],
              "FF name": [ "CPBB", "TC3", "CPNH", "CPCO" ],
              "coords": np.array([[ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00],
                                  [ 2.83557758e+00,  5.34189583e-08,  3.07968850e-07],
                                  [ 0.00000000e+00,  3.50000000e-01,  1.05300000e+00],
                                  [ 0.00000000e+00, -3.50000000e-01,  1.05300000e+00]]),
              "charge": [ 0, 0, 0.4, -0.4 ],
              "BB bonds": [],
              "SC bonds": [],
              "constraints": np.array([[ 0, 1, 0.27]]),
              "BB angles": [],
              "SC angles": [],
              "BB dihedral": [],
              "SC dihedral": [],
              "BB ref": [ 0 ],
              "BB-SC ref": [ 0, 1 ],
              "BB-SC angle": 100,
              "vsites3": { 2: [ 0.35, 0.00, 1.053 ],
                           3: [ 0.00, 0.35, 1.053 ] },
              "exclusion": [ 2, 3 ],
              "int exclusion": False,
              "vsitesn": []
              },
          "martini22p": {
              
              } 
        },
        
    # Valine
    "V" :
        { "name": "VAL",
          "AA" : {
              "beads": [ "N", "H", "CA", "CB", "CG1", "CG2", "C", "O" ],
              "coords": np.array([[ -7.071e-01, -1.225e+00,  5.000e-01 ],
                                  [ -7.071e-01, -1.225e+00,  1.500e+00 ],
                                  [  0.000e+00,  0.000e+00,  0.000e+00 ],
                                  [  1.520e+00, -3.675e-08,  5.067e-01 ],
                                  [  2.060e+00, -1.141e+00, -3.523e-01 ],
                                  [  2.052e+00,  1.312e+00, -6.497e-02 ],
                                  [ -7.165e-01,  1.241e+00,  5.067e-01 ],
                                  [ -7.165e-01,  1.241e+00,  1.737e+00 ]]),
              "BB ref": [ 0, 2, 6, 7 ],
              "SC ref": [ 3 ],
              },
          "martini22":  {
              "beads": [ "BB", "SC1" ],
              "mass": [ None, None ],
              "FF name": [ "P5", "AC2" ],
              "coords": np.array([[ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00],
                                  [ 2.83557758e+00,  5.34189583e-08,  3.07968850e-07]]),
              "charge": [ 0, 0 ],
              "BB bonds": [],
              "SC bonds": [],
              "constraints": np.array([[ 0, 1, 0.265]]),
              "BB angles": [],
              "SC angles": [],
              "BB dihedral": [],
              "SC dihedral": [],
              "BB ref": [ 0 ],
              "BB-SC ref": [ 0, 1 ],
              "BB-SC angle": 100,
              "vsites3": { },
              "exclusion": [ ],
              "int exclusion": False,
              "vsitesn": []
              },
          "marstini22": {
              "beads": [ "BB", "SC1" ,"H", "O" ],
              "mass": [ None, None, 0.00, 0.00 ],
              "FF name": [ "CPBB", "AC2", "CPNH", "CPCO" ],
              "coords": np.array([[ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00],
                                  [ 2.83557758e+00,  5.34189583e-08,  3.07968850e-07],
                                  [ 0.00000000e+00,  3.50000000e-01,  1.05300000e+00],
                                  [ 0.00000000e+00, -3.50000000e-01,  1.05300000e+00]]),
              "charge": [ 0, 0, 0.4, -0.4 ],
              "BB bonds": [],
              "SC bonds": [],
              "constraints": np.array([[ 0, 1, 0.265]]),
              "BB angles": [],
              "SC angles": [],
              "BB dihedral": [],
              "SC dihedral": [],
              "BB ref": [ 0 ],
              "BB-SC ref": [ 0, 1 ],
              "BB-SC angle": 100,
              "vsites3": { 2: [ 0.35, 0.00, 1.053 ],
                           3: [ 0.00, 0.35, 1.053 ] },
              "exclusion": [ 2, 3 ],
              "int exclusion": False,
              "vsitesn": []
              },
          "martini3": {
              "beads": [ "BB", "SC1" ],
              "mass": [ None, None ],
              "FF name": [ "P2", "SC3"],
              "coords": np.array([[ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00],
                                  [ 2.83557758e+00,  5.34189583e-08,  3.07968850e-07]]),
              "charge": [ 0, 0 ],
              "BB bonds": [],
              "SC bonds": [],
              "constraints": np.array([[ 0, 1, 0.292]]),
              "BB angles": [],
              "SC angles": [],
              "BB dihedral": [],
              "SC dihedral": [],
              "BB ref": [ 0 ],
              "BB-SC ref": [ 0, 1 ],
              "BB-SC angle": 100,
              "vsites3": {},
              "exclusion": [ ],
              "int exclusion": False,
              "vsitesn": []
              },
          "marstini3": {
              "beads": [ "BB", "SC1" ,"H", "O" ],
              "mass": [ None, None, 0.00, 0.00 ],
              "FF name": [ "CPBB", "SC3", "CPNH", "CPCO" ],
              "coords": np.array([[ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00],
                                  [ 2.83557758e+00,  5.34189583e-08,  3.07968850e-07],
                                  [ 0.00000000e+00,  3.50000000e-01,  1.05300000e+00],
                                  [ 0.00000000e+00, -3.50000000e-01,  1.05300000e+00]]),
              "charge": [ 0, 0, 0.4, -0.4 ],
              "BB bonds": [],
              "SC bonds": [],
              "constraints": np.array([[ 0, 1, 0.292]]),
              "BB angles": [],
              "SC angles": [],
              "BB dihedral": [],
              "SC dihedral": [],
              "BB ref": [ 0 ],
              "BB-SC ref": [ 0, 1 ],
              "BB-SC angle": 100,
              "vsites3": { 2: [ 0.35, 0.00, 1.053 ],
                           3: [ 0.00, 0.35, 1.053 ] },
              "exclusion": [ 2, 3 ],
              "int exclusion": False,
              "vsitesn": []
              },
          "martini22p": {
              
              } 
        },
        
    # Isoleucine
    "I" :
        { "name": "ILE",
          "AA" : {
              "beads": [ "N", "H", "CA", "CB", "CG1", "CG2", "CD", "C", "O" ],
              "coords": np.array([[ -7.071e-01, -1.225e+00,  5.000e-01 ],
                                  [ -7.071e-01, -1.225e+00,  1.500e+00 ],
                                  [  0.000e+00,  0.000e+00,  0.000e+00 ],
                                  [  1.520e+00, -1.805e-08,  5.067e-01 ],
                                  [  2.060e+00,  1.055e+00, -4.575e-01 ],
                                  [  2.053e+00, -1.374e+00,  1.048e-01 ],
                                  [  3.556e+00,  1.290e+00, -3.396e-01 ],
                                  [ -7.165e-01,  1.241e+00,  5.067e-01 ],
                                  [ -7.165e-01,  1.241e+00,  1.737e+00 ]]),
              "BB ref": [ 0, 2, 7, 8 ],
              "SC ref": [ 3 ],
              },
          "martini22":  {
              "beads": [ "BB", "SC1" ],
              "mass": [ None, None ],
              "FF name": [ "P5", "AC1" ],
              "coords": np.array([[ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00],
                                  [ 2.99153782e+00, -1.55740071e-09, -7.41370540e-09]]),
              "charge": [ 0, 0 ],
              "BB bonds": [],
              "SC bonds": [],
              "constraints": np.array([[ 0, 1, 0.31]]),
              "BB angles": [],
              "SC angles": [],
              "BB dihedral": [],
              "SC dihedral": [],
              "BB ref": [ 0 ],
              "BB-SC ref": [ 0, 1 ],
              "BB-SC angle": 100,
              "vsites3": { },
              "exclusion": [ ],
              "int exclusion": False,
              "vsitesn": []
              },
          "marstini22": {
              "beads": [ "BB", "SC1", "H", "O" ],
              "mass": [ None, None, 0.00, 0.00 ],
              "FF name": [ "CPBB", "AC1", "CPNH", "CPCO" ],
              "coords": np.array([[ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00],
                                  [ 2.99153782e+00, -1.55740071e-09, -7.41370540e-09],
                                  [ 0.00000000e+00,  3.50000000e-01,  1.05300000e+00],
                                  [ 0.00000000e+00, -3.50000000e-01,  1.05300000e+00]]),
              "charge": [ 0, 0, 0.4, -0.4 ],
              "BB bonds": [],
              "SC bonds": [],
              "constraints": np.array([[ 0, 1, 0.31 ]]),
              "BB angles": [],
              "SC angles": [],
              "BB dihedral": [],
              "SC dihedral": [],
              "BB ref": [ 0 ],
              "BB-SC ref": [ 0, 1 ],
              "BB-SC angle": 100,
              "vsites3": { 2: [ 0.35, 0.00, 1.053 ],
                           3: [ 0.00, 0.35, 1.053 ] },
              "exclusion": [ 2, 3 ],
              "int exclusion": False,
              "vsitesn": []
              },
          "martini3": {
              "beads": [ "BB", "SC1"],
              "mass": [ None, None ],
              "FF name": [ "P2", "SC2"],
              "coords": np.array([[ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00],
                                  [ 2.99153782e+00, -1.55740071e-09, -7.41370540e-09]]),
              "charge": [ 0, 0],
              "BB bonds": [],
              "SC bonds": [],
              "constraints": np.array([[ 0, 1, 0.341 ]]),
              "BB angles": [],
              "SC angles": [],
              "BB dihedral": [],
              "SC dihedral": [],
              "BB ref": [ 0 ],
              "BB-SC ref": [ 0, 1 ],
              "BB-SC angle": 100,
              "vsites3": {},
              "exclusion": [ ],
              "int exclusion": False,
              "vsitesn": []
              },
          "marstini3": {
              "beads": [ "BB", "SC1", "H", "O" ],
              "mass": [ None, None, 0.00, 0.00 ],
              "FF name": [ "CPBB", "SC2", "CPNH", "CPCO" ],
              "coords": np.array([[ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00],
                                  [ 2.99153782e+00, -1.55740071e-09, -7.41370540e-09],
                                  [ 0.00000000e+00,  3.50000000e-01,  1.05300000e+00],
                                  [ 0.00000000e+00, -3.50000000e-01,  1.05300000e+00]]),
              "charge": [ 0, 0, 0.4, -0.4 ],
              "BB bonds": [],
              "SC bonds": [],
              "constraints": np.array([[ 0, 1, 0.341 ]]),
              "BB angles": [],
              "SC angles": [],
              "BB dihedral": [],
              "SC dihedral": [],
              "BB ref": [ 0 ],
              "BB-SC ref": [ 0, 1 ],
              "BB-SC angle": 100,
              "vsites3": { 2: [ 0.35, 0.00, 1.053 ],
                           3: [ 0.00, 0.35, 1.053 ] },
              "exclusion": [ 2, 3 ],
              "int exclusion": False,
              "vsitesn": []
              },
          "martini22p": {
              
              } 
        },
        
    # Leucine
    "L" :
        { "name": "LEU",
          "AA" : {
              "beads": [ "N", "H", "CA", "CB", "CG", "CD1", "CD2", "C", "O" ],
              "coords": np.array([[ -7.070e-01, -1.225e+00,  5.000e-01 ],
                                  [ -7.070e-01, -1.225e+00,  1.500e+00 ],
                                  [  0.000e+00,  0.000e+00,  0.000e+00 ],
                                  [  1.520e+00,  0.000e+00,  5.070e-01 ],
                                  [  2.193e+00,  1.374e+00,  5.070e-01 ],
                                  [  3.704e+00,  1.221e+00,  6.320e-01 ],
                                  [  1.922e+00,  2.103e+00,  1.819e+00 ],
                                  [ -7.160e-01,  1.241e+00,  5.070e-01 ],
                                  [ -7.160e-01,  1.241e+00,  1.737e+00 ]]),
              "BB ref": [ 0, 2, 7, 8 ],
              "SC ref": [ 3 ],
              },
          "martini22":  {
              "beads": [ "BB", "SC1" ],
              "mass": [ None, None ],
              "FF name": [ "P5", "AC1" ],
              "coords": np.array([[ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00],
                                  [ 3.49905756e+00,  1.93707215e-08, -2.48815373e-08]]),
              "charge": [ 0, 0 ],
              "BB bonds": [],
              "SC bonds": np.array([[ 0, 1, 0.33, 7500 ]]),
              "constraints": [],
              "BB angles": [],
              "SC angles": [],
              "BB dihedral": [],
              "SC dihedral": [],
              "BB ref": [ 0 ],
              "BB-SC ref": [ 0, 1 ],
              "BB-SC angle": 100,
              "vsites3": { },
              "exclusion": [ ],
              "int exclusion": False,
              "vsitesn": []
              },
          "marstini22":  {
              "beads": [ "BB", "SC1", "H", "O" ],
              "mass": [ None, None, 0.00, 0.00 ],
              "FF name": [ "CPBB", "AC1", "CPNH", "CPCO" ],
              "coords": np.array([[ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00],
                                  [ 3.49905756e+00,  1.93707215e-08, -2.48815373e-08],
                                  [ 0.00000000e+00,  3.50000000e-01,  1.05300000e+00],
                                  [ 0.00000000e+00, -3.50000000e-01,  1.05300000e+00]]),
              "charge": [ 0, 0, 0.4, -0.4 ],
              "BB bonds": [],
              "SC bonds": np.array([[ 0, 1, 0.33, 7500 ]]),
              "constraints": [],
              "BB angles": [],
              "SC angles": [],
              "BB dihedral": [],
              "SC dihedral": [],
              "BB ref": [ 0 ],
              "BB-SC ref": [ 0, 1 ],
              "BB-SC angle": 100,
              "vsites3": { 2: [ 0.35, 0.00, 1.053 ],
                            3: [ 0.00, 0.35, 1.053 ] },
              "exclusion": [ 2, 3 ],
              "int exclusion": False,
              "vsitesn": []
              },
          "martini3": {
              "beads": [ "BB", "SC1" ],
              "mass": [ None, None ],
              "FF name": [ "P2", "SC2" ],
              "coords": np.array([[ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00],
                                  [ 3.49905756e+00,  1.93707215e-08, -2.48815373e-08]]),
              "charge": [ 0, 0 ],
              "BB bonds": [],
              "SC bonds": np.array([[ 0, 1, 0.363, 7500 ]]),
              "constraints": [],
              "BB angles": [],
              "SC angles": [],
              "BB dihedral": [],
              "SC dihedral": [],
              "BB ref": [ 0 ],
              "BB-SC ref": [ 0, 1 ],
              "BB-SC angle": 100,
              "vsites3": { },
              "exclusion": [ ],
              "int exclusion": False,
              "vsitesn": []
              },
          "marstini3": {
              "beads": [ "BB", "SC1", "H", "O" ],
              "mass": [ None, None, 0.00, 0.00 ],
              "FF name": [ "CPBB", "SC2", "CPNH", "CPCO" ],
              "coords": np.array([[ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00],
                                  [ 3.49905756e+00,  1.93707215e-08, -2.48815373e-08],
                                  [ 0.00000000e+00,  3.50000000e-01,  1.05300000e+00],
                                  [ 0.00000000e+00, -3.50000000e-01,  1.05300000e+00]]),
              "charge": [ 0, 0, 0.4, -0.4 ],
              "BB bonds": [],
              "SC bonds": np.array([[ 0, 1, 0.363, 7500 ]]),
              "constraints": [],
              "BB angles": [],
              "SC angles": [],
              "BB dihedral": [],
              "SC dihedral": [],
              "BB ref": [ 0 ],
              "BB-SC ref": [ 0, 1 ],
              "BB-SC angle": 100,
              "vsites3": { 2: [ 0.35, 0.00, 1.053 ],
                            3: [ 0.00, 0.35, 1.053 ] },
              "exclusion": [ 2, 3 ],
              "int exclusion": False,
              "vsitesn": []
              },
          "martini22p": {
              
              } 
        },
        
    # Methionine
    "M" :
        { "name": "MET",
          "AA" : {
              "beads": [ "N", "H", "CA", "CB", "CG", "SD", "CE", "C", "O" ],
              "coords": np.array([[ -7.070e-01, -1.225e+00,  5.000e-01 ],
                                  [ -7.070e-01, -1.225e+00,  1.500e+00 ],
                                  [  0.000e+00,  0.000e+00,  0.000e+00 ],
                                  [  1.520e+00,  0.000e+00,  5.070e-01 ],
                                  [  2.131e+00,  1.392e+00,  5.070e-01 ],
                                  [  3.941e+00,  1.360e+00,  4.960e-01 ],
                                  [  4.239e+00, -4.040e-01,  4.900e-01 ],
                                  [ -7.160e-01,  1.241e+00,  5.070e-01 ],
                                  [ -7.160e-01,  1.241e+00,  1.737e+00 ]]),
              "BB ref": [ 0, 2, 7, 8 ],
              "SC ref": [ 3 ],
              },
          "martini22":  {
              "beads": [ "BB", "SC1" ],
              "mass": [ None, None ],
              "FF name": [ "P5", "C5" ],
              "coords": np.array([[ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00],
                                  [ 4.31197328e+00,  1.13892811e-07, -2.87471216e-07]]),
              "charge": [ 0, 0 ],
              "BB bonds": [],
              "SC bonds": np.array([[ 0, 1, 0.40, 2500 ]]),
              "constraints": [],
              "BB angles": [],
              "SC angles": [],
              "BB dihedral": [],
              "SC dihedral": [],
              "BB ref": [ 0 ],
              "BB-SC ref": [ 0, 1 ],
              "BB-SC angle": 100,
              "vsites3": { },
              "exclusion": [ ],
              "int exclusion": False,
              "vsitesn": []
              },
          "marstini22": {
              "beads": [ "BB", "SC1", "H", "O" ],
              "mass": [ None, None, 0.00, 0.00 ],
              "FF name": [ "CPBB", "C5", "CPNH", "CPCO" ],
              "coords": np.array([[ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00],
                                  [ 4.31197328e+00,  1.13892811e-07, -2.87471216e-07],
                                  [ 0.00000000e+00,  3.50000000e-01,  1.05300000e+00],
                                  [ 0.00000000e+00, -3.50000000e-01,  1.05300000e+00]]),
              "charge": [ 0, 0, 0.4, -0.4 ],
              "BB bonds": [],
              "SC bonds": np.array([[ 0, 1, 0.40, 2500 ]]),
              "constraints": [],
              "BB angles": [],
              "SC angles": [],
              "BB dihedral": [],
              "SC dihedral": [],
              "BB ref": [ 0 ],
              "BB-SC ref": [ 0, 1 ],
              "BB-SC angle": 100,
              "vsites3": { 2: [ 0.35, 0.00, 1.053 ],
                            3: [ 0.00, 0.35, 1.053 ] },
              "exclusion": [ 2, 3 ],
              "int exclusion": False,
              "vsitesn": []
              },
          "martini3": {
              "beads": [ "BB", "SC1" ],
              "mass": [ None, None ],
              "FF name": [ "P2", "C6" ],
              "coords": np.array([[ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00],
                                  [ 4.31197328e+00,  1.13892811e-07, -2.87471216e-07]]),
              "charge": [ 0, 0 ],
              "BB bonds": [],
              "SC bonds": np.array([[ 0, 1, 0.40, 2500 ]]),
              "constraints": [],
              "BB angles": [],
              "SC angles": [],
              "BB dihedral": [],
              "SC dihedral": [],
              "BB ref": [ 0 ],
              "BB-SC ref": [ 0, 1 ],
              "BB-SC angle": 100,
              "vsites3": {},
              "exclusion": [ ],
              "int exclusion": False,
              "vsitesn": []
              },
          "marstini3": {
              "beads": [ "BB", "SC1", "H", "O" ],
              "mass": [ None, None, 0.00, 0.00 ],
              "FF name": [ "CPBB", "C6", "CPNH", "CPCO" ],
              "coords": np.array([[ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00],
                                  [ 4.31197328e+00,  1.13892811e-07, -2.87471216e-07],
                                  [ 0.00000000e+00,  3.50000000e-01,  1.05300000e+00],
                                  [ 0.00000000e+00, -3.50000000e-01,  1.05300000e+00]]),
              "charge": [ 0, 0, 0.4, -0.4 ],
              "BB bonds": [],
              "SC bonds": np.array([[ 0, 1, 0.40, 2500 ]]),
              "constraints": [],
              "BB angles": [],
              "SC angles": [],
              "BB dihedral": [],
              "SC dihedral": [],
              "BB ref": [ 0 ],
              "BB-SC ref": [ 0, 1 ],
              "BB-SC angle": 100,
              "vsites3": { 2: [ 0.35, 0.00, 1.053 ],
                            3: [ 0.00, 0.35, 1.053 ] },
              "exclusion": [ 2, 3 ],
              "int exclusion": False,
              "vsitesn": []
              },
          "martini22p": {
              
              } 
        },
        
    # Phenylalanine
    "F" :
        { "name": "PHE",
          "AA" : {
              "beads": [ "N", "H", "CA", "CB", "CG", "CD1", "HD1", "CD2", "HD2", "CE1", "HE1", "CE2", "HE2", "CZ", "HZ", "C", "O" ],
              "coords": np.array([[ -7.070e-01, -1.225e+00,  5.000e-01 ],
                                  [ -7.070e-01, -1.225e+00,  1.500e+00 ],
                                  [  0.000e+00,  0.000e+00,  0.000e+00 ],
                                  [  1.521e+00, -0.000e+00,  5.070e-01 ],
                                  [  2.127e+00,  1.372e+00,  5.070e-01 ],
                                  [  1.306e+00,  2.494e+00,  5.070e-01 ],
                                  [  3.120e-01,  2.387e+00,  5.040e-01 ],
                                  [  3.508e+00,  1.521e+00,  5.090e-01 ],
                                  [  4.099e+00,  7.130e-01,  5.080e-01 ],
                                  [  1.868e+00,  3.765e+00,  5.080e-01 ],
                                  [  1.278e+00,  4.572e+00,  5.090e-01 ],
                                  [  4.070e+00,  2.793e+00,  5.100e-01 ],
                                  [  5.065e+00,  2.901e+00,  5.120e-01 ],
                                  [  3.250e+00,  3.915e+00,  5.090e-01 ],
                                  [  3.654e+00,  4.830e+00,  5.080e-01 ],
                                  [ -7.160e-01,  1.241e+00,  5.070e-01 ],
                                  [ -7.160e-01,  1.241e+00,  1.737e+00 ]]),
              "BB ref": [ 0, 2, 15, 16 ],
              "SC ref": [ 3 ],
              },
          "martini22":  {
              "beads": [ "BB", "SC1", "SC2", "SC3" ],
              "mass": [ None, None, None, None ],
              "FF name": [ "P5", "SC5", "SC5", "SC5" ],
              "coords": np.array([[ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00],
                                  [ 3.24336542e+00, -1.21474470e-07, -1.36387067e-07],
                                  [ 4.65658902e+00,  5.65779030e-01,  1.73380690e+00],
                                  [ 5.55109274e+00, -1.10511949e+00,  8.87013314e-01]]),
              "charge": [ 0, 0, 0, 0 ],
              "BB bonds": [],
              "SC bonds": np.array([[ 0, 1, 0.31, 7500 ]]),
              "constraints": np.array([[ 1, 2, 0.27 ],
                                       [ 1, 3, 0.27 ],
                                       [ 2, 3, 0.27 ]]),
              "BB angles": [],
              "SC angles": np.array([[ 0, 1, 2, 150, 50 ],
                                     [ 0, 1, 3, 150, 50 ]]),
              "BB dihedral": [],
              "SC dihedral": np.array([[ 0, 2, 3, 1, 0, 50 ]]),
              "BB ref": [ 0 ],
              "BB-SC ref": [ 0, 1 ],
              "BB-SC angle": 100,
              "vsites3": { },
              "exclusion": [ ],
              "int exclusion": False,
              "vsitesn": []
              },
          "marstini22": {
              "beads": [ "BB", "SC1", "SC2", "SC3", "H", "O" ],
              "mass": [ None, None, None, None, 0.00, 0.00 ],
              "FF name": [ "CPBB", "SC5", "SC5", "SC5", "CPNH", "CPCO" ],
              "coords": np.array([[ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00],
                                  [ 3.24336542e+00, -1.21474470e-07, -1.36387067e-07],
                                  [ 4.65658902e+00,  5.65779030e-01,  1.73380690e+00],
                                  [ 5.55109274e+00, -1.10511949e+00,  8.87013314e-01],
                                  [ 0.00000000e+00,  3.50000000e-01,  1.05300000e+00],
                                  [ 0.00000000e+00, -3.50000000e-01,  1.05300000e+00]]),
              "charge": [ 0, 0, 0, 0, 0.4, -0.4 ],
              "BB bonds": [],
              "SC bonds": np.array([[ 0, 1, 0.31, 7500 ]]),
              "constraints": np.array([[ 1, 2, 0.27 ],
                                       [ 1, 3, 0.27 ],
                                       [ 2, 3, 0.27 ]]),
              "BB angles": [],
              "SC angles": np.array([[ 0, 1, 2, 150, 50 ],
                                     [ 0, 1, 3, 150, 50 ]]),
              "BB dihedral": [],
              "SC dihedral": np.array([[ 0, 2, 3, 1, 0, 50 ]]),
              "BB ref": [ 0 ],
              "BB-SC ref": [ 0, 1 ],
              "BB-SC angle": 100,
              "vsites3": { 4: [ 0.35, 0.00, 1.053 ],
                           5: [ 0.00, 0.35, 1.053 ] },
              "exclusion": [ 4, 5 ],
              "int exclusion": False,
              "vsitesn": []
              },
          "martini3": {
              "beads": [ "BB", "SC1", "SC2", "SC3"],
              "mass": [ None, None, None, None ],
              "FF name": [ "P2", "SC4", "TC5", "TC5"],
              "coords": np.array([[ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00],
                                  [ 3.24336542e+00, -1.21474470e-07, -1.36387067e-07],
                                  [ 4.65658902e+00,  5.65779030e-01,  1.73380690e+00],
                                  [ 5.55109274e+00, -1.10511949e+00,  8.87013314e-01]]),
              "charge": [ 0, 0, 0, 0 ],
              "BB bonds": [],
              "SC bonds": np.array([[ 0, 1, 0.325, 7500 ]]),
              "constraints": np.array([[ 1, 2, 0.34 ],
                                       [ 1, 3, 0.34 ],
                                       [ 2, 3, 0.29 ]]),
              "BB angles": [],
              "SC angles": np.array([[ 0, 1, 2, 120, 50 ],
                                     [ 0, 1, 3, 120, 50 ]]),
              "BB dihedral": [],
              "SC dihedral": [],
              "BB ref": [ 0 ],
              "BB-SC ref": [ 0, 1 ],
              "BB-SC angle": 100,
              "vsites3": { },
              "exclusion": [ ],
              "int exclusion": True,
              "vsitesn": []
              },
          "marstini3": {
              "beads": [ "BB", "SC1", "SC2", "SC3", "H", "O" ],
              "mass": [ None, None, None, None, 0.00, 0.00 ],
              "FF name": [ "CPBB", "SC4", "TC5", "TC5", "CPNH", "CPCO" ],
              "coords": np.array([[ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00],
                                  [ 3.24336542e+00, -1.21474470e-07, -1.36387067e-07],
                                  [ 4.65658902e+00,  5.65779030e-01,  1.73380690e+00],
                                  [ 5.55109274e+00, -1.10511949e+00,  8.87013314e-01],
                                  [ 0.00000000e+00,  3.50000000e-01,  1.05300000e+00],
                                  [ 0.00000000e+00, -3.50000000e-01,  1.05300000e+00]]),
              "charge": [ 0, 0, 0, 0, 0.4, -0.4 ],
              "BB bonds": [],
              "SC bonds": np.array([[ 0, 1, 0.325, 7500 ]]),
              "constraints": np.array([[ 1, 2, 0.34 ],
                                       [ 1, 3, 0.34 ],
                                       [ 2, 3, 0.29 ]]),
              "BB angles": [],
              "SC angles": np.array([[ 0, 1, 2, 120, 50 ],
                                     [ 0, 1, 3, 120, 50 ]]),
              "BB dihedral": [],
              "SC dihedral": [],
              "BB ref": [ 0 ],
              "BB-SC ref": [ 0, 1 ],
              "BB-SC angle": 100,
              "vsites3": { 4: [ 0.35, 0.00, 1.053 ],
                           5: [ 0.00, 0.35, 1.053 ] },
              "exclusion": [ 4, 5 ],
              "int exclusion": True,
              "vsitesn": []
              },
          "martini22p": {
              
              } 
        },
        
    # Tyrosine
    "Y" :
        { "name": "TYR",
          "AA" : {
              "beads": [ "N", "H", "CA", "CB", "CG", "CD1", "HD1", "CD2", "HD2", "CE1", "HE1", "CE2", "HE2", "CZ", "OH", "HH", "C", "O" ],
              "coords": np.array([[ -7.070e-01, -1.225e+00,  5.000e-01 ],
                                  [ -7.070e-01, -1.225e+00,  1.500e+00 ],
                                  [  0.000e+00,  0.000e+00,  0.000e+00 ],
                                  [  1.520e+00,  0.000e+00,  5.070e-01 ],
                                  [  2.130e+00,  1.381e+00,  5.070e-01 ],
                                  [  1.328e+00,  2.517e+00,  5.070e-01 ],
                                  [  3.330e-01,  2.426e+00,  5.080e-01 ],
                                  [  3.509e+00,  1.551e+00,  5.070e-01 ],
                                  [  4.111e+00,  7.530e-01,  5.070e-01 ],
                                  [  1.910e+00,  3.779e+00,  5.070e-01 ],
                                  [  1.334e+00,  4.596e+00,  5.070e-01 ],
                                  [  4.051e+00,  2.830e+00,  5.070e-01 ],
                                  [  5.046e+00,  2.937e+00,  5.070e-01 ],
                                  [  3.295e+00,  3.905e+00,  5.060e-01 ],
                                  [  3.872e+00,  5.169e+00,  5.060e-01 ],
                                  [  4.868e+00,  5.081e+00,  5.060e-01 ],
                                  [ -7.160e-01,  1.241e+00,  5.070e-01 ],
                                  [ -7.160e-01,  1.241e+00,  1.737e+00 ]]),
              "BB ref": [ 0, 2, 16, 17 ],
              "SC ref": [ 3 ],
              },
          "martini22":  {
              "beads": [ "BB", "SC1", "SC2", "SC3" ],
              "mass": [ None, None, None, None ],
              "FF name": [ "P5", "SC4", "SC4", "SP1" ],
              "coords": np.array([[ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00],
                                  [ 3.24919957e+00, -5.42288713e-08,  5.64100793e-08],
                                  [ 4.67241023e+00, -1.69478236e+00,  5.96975767e-01],
                                  [ 6.15499141e+00, -1.37744538e+00, -1.10109824e+00]]),
              "charge": [ 0, 0, 0, 0 ],
              "BB bonds": [],
              "SC bonds": np.array([[ 0, 1, 0.32, 5000 ]]),
              "constraints": np.array([[ 1, 2, 0.27 ],
                                       [ 1, 3, 0.27 ],
                                       [ 2, 3, 0.27 ]]),
              "BB angles": [],
              "SC angles": np.array([[ 0, 1, 2, 150, 50 ],
                                     [ 0, 1, 3, 150, 50 ]]),
              "BB dihedral": [],
              "SC dihedral": np.array([[ 0, 2, 3, 1, 0, 50 ]]),
              "BB ref": [ 0 ],
              "BB-SC ref": [ 0, 1 ],
              "BB-SC angle": 100,
              "vsites3": { },
              "exclusion": [ ],
              "int exclusion": False,
              "vsitesn": []
              },
          "marstini22": {
              "beads": [ "BB", "SC1", "SC2", "SC3", "H", "O" ],
              "mass": [ None, None, None, None, 0.00, 0.00 ],
              "FF name": [ "CPBB", "SC4", "SC4", "SP1", "CPNH", "CPCO" ],
              "coords": np.array([[ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00],
                                  [ 3.24919957e+00, -5.42288713e-08,  5.64100793e-08],
                                  [ 4.67241023e+00, -1.69478236e+00,  5.96975767e-01],
                                  [ 6.15499141e+00, -1.37744538e+00, -1.10109824e+00],
                                  [ 0.00000000e+00,  3.50000000e-01,  1.05300000e+00],
                                  [ 0.00000000e+00, -3.50000000e-01,  1.05300000e+00]]),
              "charge": [ 0, 0, 0, 0, 0.4, -0.4 ],
              "BB bonds": [],
              "SC bonds": np.array([[ 0, 1, 0.32, 5000 ]]),
              "constraints": np.array([[ 1, 2, 0.27 ],
                                       [ 1, 3, 0.27 ],
                                       [ 2, 3, 0.27 ]]),
              "BB angles": [],
              "SC angles": np.array([[ 0, 1, 2, 150, 50 ],
                                     [ 0, 1, 3, 150, 50 ]]),
              "BB dihedral": [],
              "SC dihedral": np.array([[ 0, 2, 3, 1, 0, 50 ]]),
              "BB ref": [ 0 ],
              "BB-SC ref": [ 0, 1 ],
              "BB-SC angle": 100,
              "vsites3": { 4: [ 0.35, 0.00, 1.053 ],
                           5: [ 0.00, 0.35, 1.053 ] },
              "exclusion": [ 4, 5 ],
              "int exclusion": False,
              "vsitesn": []
              },
          "martini3": {
              "beads": [ "BB", "SC1", "SC2", "SC3", "SC4" ],
              "mass": [ None, None, None, None, None ],
              "FF name": [ "P2", "TC4", "TC5", "TC5", "TN6"],
              "coords": np.array([[ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00],
                                  [ 3.36413405e+00, -3.30197799e-07, -3.98250631e-07],
                                  [ 4.87222078e+00, -1.44293460e+00, -3.49448453e-02],
                                  [ 5.02179871e+00,  1.15136211e+00,  2.13923300e+00],
                                  [ 6.52114639e+00, -2.22384403e-01,  2.16703598e+00]]),
              "charge": [ 0, 0, 0, 0, 0 ],
              "BB bonds": [],
              "SC bonds": np.array([[ 0, 1, 0.325, 5000 ]]),
              "constraints": np.array([[ 1, 2, 0.30 ],
                                       [ 1, 3, 0.30 ],
                                       [ 2, 3, 0.3 ],
                                       [ 2, 4, 0.285],
                                       [ 3, 4, 0.285]]),
              "BB angles": [],
              "SC angles": np.array([[ 0, 1, 2, 120, 60 ],
                                     [ 0, 1, 3, 120, 60 ]]),
              "BB dihedral": [],
              "SC dihedral": np.array([[ 4, 2, 3, 1, 180, 50 ]]),
              "BB ref": [ 0 ],
              "BB-SC ref": [ 0, 1 ],
              "BB-SC angle": 100,
              "vsites3": { },
              "exclusion": [ ],
              "int exclusion": True,
              "vsitesn": []
              },
          "marstini3": {
              "beads": [ "BB", "SC1", "SC2", "SC3", "SC4", "H", "O" ],
              "mass": [ None, None, None, None, None, 0.00, 0.00 ],
              "FF name": [ "CPBB", "TC4", "TC5", "TC5", "TN6", "CPNH", "CPCO" ],
              "coords": np.array([[ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00],
                                  [ 3.36413405e+00, -3.30197799e-07, -3.98250631e-07],
                                  [ 4.87222078e+00, -1.44293460e+00, -3.49448453e-02],
                                  [ 5.02179871e+00,  1.15136211e+00,  2.13923300e+00],
                                  [ 6.52114639e+00, -2.22384403e-01,  2.16703598e+00],
                                  [ 0.00000000e+00,  3.50000000e-01,  1.05300000e+00],
                                  [ 0.00000000e+00, -3.50000000e-01,  1.05300000e+00]]),
              "charge": [ 0, 0, 0, 0, 0, 0.4, -0.4 ],
              "BB bonds": [],
              "SC bonds": np.array([[ 0, 1, 0.325, 5000 ]]),
              "constraints": np.array([[ 1, 2, 0.30 ],
                                       [ 1, 3, 0.30 ],
                                       [ 2, 3, 0.3 ],
                                       [ 2, 4, 0.285],
                                       [ 3, 4, 0.285]]),
              "BB angles": [],
              "SC angles": np.array([[ 0, 1, 2, 120, 60 ],
                                     [ 0, 1, 3, 120, 60 ]]),
              "BB dihedral": [],
              "SC dihedral": np.array([[ 4, 2, 3, 1, 180, 50 ]]),
              "BB ref": [ 0 ],
              "BB-SC ref": [ 0, 1 ],
              "BB-SC angle": 100,
              "vsites3": { 5: [ 0.35, 0.00, 1.053 ],
                           6: [ 0.00, 0.35, 1.053 ] },
              "exclusion": [ 5, 6 ],
              "int exclusion": True,
              "vsitesn": []
              },
          "martini22p": {
              
              } 
        },
        
    # Tryptophan
    "W" :
        { "name": "TRP",
          "AA": {
              "beads": [ "N", "H", "CA", "CB", "CG", "CD1", "HD1", "CD2", "NE1", "HE1", "CE2", "CE3", "HE3", "CZ2", "HZ2", "CZ3", "HZ3", "CH2", "HH2", "C", "O" ],
              "coords": np.array([[ -7.070e-01, -1.225e+00,  5.000e-01 ],
                                  [ -7.070e-01, -1.225e+00,  1.500e+00 ],
                                  [  0.000e+00,  0.000e+00,  0.000e+00 ],
                                  [  1.520e+00, -0.000e+00,  5.070e-01 ],
                                  [  2.132e+00,  1.369e+00,  5.070e-01 ],
                                  [  1.470e+00,  2.569e+00,  5.070e-01 ],
                                  [  4.780e-01,  2.689e+00,  5.060e-01 ],
                                  [  3.527e+00,  1.682e+00,  5.060e-01 ],
                                  [  2.404e+00,  3.584e+00,  5.070e-01 ],
                                  [  2.185e+00,  4.559e+00,  5.080e-01 ],
                                  [  3.671e+00,  3.075e+00,  5.070e-01 ],
                                  [  4.694e+00,  9.080e-01,  5.070e-01 ],
                                  [  4.631e+00, -9.000e-02,  5.060e-01 ],
                                  [  4.949e+00,  3.646e+00,  5.070e-01 ],
                                  [  5.051e+00,  4.641e+00,  5.060e-01 ],
                                  [  5.947e+00,  1.530e+00,  5.070e-01 ],
                                  [  6.768e+00,  9.590e-01,  5.070e-01 ],
                                  [  6.083e+00,  2.825e+00,  5.070e-01 ],
                                  [  6.996e+00,  3.233e+00,  5.070e-01 ],
                                  [ -7.160e-01,  1.241e+00,  5.070e-01 ],
                                  [ -7.160e-01,  1.241e+00,  1.737e+00 ]]),
              "BB ref": [ 0, 2, 19, 20 ],
              "SC ref": [ 3 ],
              },
          "martini22":  {
              "beads": [ "BB", "SC1", "SC2", "SC3", "SC4" ],
              "mass": [ None, None, None, None, None ],
              "FF name": [ "P5", "SC4", "SNd", "SC5", "SC5" ],
              "coords": np.array([[ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00],
                                  [ 3.36413405e+00, -3.30197799e-07, -3.98250631e-07],
                                  [ 4.87222078e+00, -1.44293460e+00, -3.49448453e-02],
                                  [ 5.02179871e+00,  1.15136211e+00,  2.13923300e+00],
                                  [ 6.52114639e+00, -2.22384403e-01,  2.16703598e+00]]),
              "charge": [ 0, 0, 0, 0, 0 ],
              "BB bonds": [],
              "SC bonds": np.array([[ 0, 1, 0.30, 5000 ]]),
              "constraints": np.array([[ 1, 2, 0.27 ],
                                       [ 1, 3, 0.27 ],
                                       [ 2, 3, 0.27 ],
                                       [ 2, 4, 0.27 ],
                                       [ 3, 4, 0.27 ]]),
              "BB angles": [],
              "SC angles": np.array([[ 0, 1, 2, 210, 50 ],
                                     [ 0, 1, 3, 90, 50 ]]),
              "BB dihedral": [],
              "SC dihedral": np.array([[ 0, 2, 3, 1, 0, 50 ],
                                       [ 1, 2, 4, 3, 0, 200 ]]),
              "BB ref": [ 0 ],
              "BB-SC ref": [ 0, 1 ],
              "BB-SC angle": 100,
              "vsites3": { },
              "exclusion": [ ],
              "int exclusion": False,
              "vsitesn": []
              },
          "marstini22": {
              "beads": [ "BB", "SC1", "SC2", "SC3", "SC4", "H", "O" ],
              "mass": [ None, None, None, None, None, 0.00, 0.00 ],
              "FF name": [ "CPBB", "SC4", "SNd", "SC5", "SC5", "CPNH", "CPCO" ],
              "coords": np.array([[ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00],
                                  [ 3.36413405e+00, -3.30197799e-07, -3.98250631e-07],
                                  [ 4.87222078e+00, -1.44293460e+00, -3.49448453e-02],
                                  [ 5.02179871e+00,  1.15136211e+00,  2.13923300e+00],
                                  [ 6.52114639e+00, -2.22384403e-01,  2.16703598e+00],
                                  [ 0.00000000e+00,  3.50000000e-01,  1.05300000e+00],
                                  [ 0.00000000e+00, -3.50000000e-01,  1.05300000e+00]]),
              "charge": [ 0, 0, 0, 0, 0, 0.4, -0.4 ],
              "BB bonds": [],
              "SC bonds": np.array([[ 0, 1, 0.30, 5000 ]]),
              "constraints": np.array([[ 1, 2, 0.27 ],
                                       [ 1, 3, 0.27 ],
                                       [ 2, 3, 0.27 ],
                                       [ 2, 4, 0.27 ],
                                       [ 3, 4, 0.27 ]]),
              "BB angles": [],
              "SC angles": np.array([[ 0, 1, 2, 210, 50 ],
                                     [ 0, 1, 3, 90, 50 ]]),
              "BB dihedral": [],
              "SC dihedral": np.array([[ 0, 2, 3, 1, 0, 50 ],
                                       [ 1, 2, 4, 3, 0, 200 ]]),
              "BB ref": [ 0 ],
              "BB-SC ref": [ 0, 1 ],
              "BB-SC angle": 100,
              "vsites3": { 5: [ 0.35, 0.00, 1.053 ],
                           6: [ 0.00, 0.35, 1.053 ] },
              "exclusion": [ 5, 6 ],
              "int exclusion": False,
              "vsitesn": []
              },
          "martini3": {
              "beads": [ "BB", "SC1", "SC2", "SC3", "SC4", "SC5"],
              "mass": [ None, 36.0, 36.0, 0.00, 36.0, 36.0 ],
              "FF name": [ "P2", "TC4", "TN6d", "TC5", "TC5", "TC5" ],
              "coords": np.array([[ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00],
                                  [ 1.90585845e+00,  5.02697766e-01, -1.87914821e+00],
                                  [ 3.24316111e+00, -1.35305546e+00,  5.59683825e-01],
                                  [ 4.49272713e+00,  2.25802114e-01, -8.05304903e-01],
                                  [ 5.74466053e+00,  1.75364395e+00, -1.90801145e+00],
                                  [ 7.09622459e+00, -1.87894658e-10,  6.33205622e-08]]),
              "charge": [ 0, 0, 0, 0, 0, 0],
              "BB bonds": [],
              "SC bonds": np.array([[ 0, 1, 0.315, 5000 ]]),
              "constraints": np.array([[ 1, 2, 0.335],
                                       [ 2, 5, 0.412],
                                       [ 4, 5, 0.293],
                                       [ 1, 4, 0.404 ],
                                       [ 2, 4, 0.470 ]]),
              "BB angles": [],
              "SC angles": np.array([[ 0, 1, 2, 120, 60 ],
                                     [ 0, 1, 4, 130, 60 ]]),
              "BB dihedral": [],
              "SC dihedral": np.array([[ 5, 4, 2, 1, 180, 100]]),
              "BB ref": [ 0 ],
              "BB-SC ref": [ 0, 1 ],
              "BB-SC angle": 100,
              "vsites3": { },
              "exclusion": [ ],
              "int exclusion": True,
              "vsitesn": [3]
              },
          "marstini3": {
              "beads": [ "BB", "SC1", "SC2", "SC3", "SC4", "SC5", "H", "O" ],
              "mass": [ None, 36.0, 36.0, 0.00, 36.0, 36.0, 0.00, 0.00 ],
              "FF name": [ "CPBB", "TC4", "TN6d", "TC5", "TC5", "TC5", "CPNH", "CPCO" ],
              "coords": np.array([[ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00],
                                  [ 1.90585845e+00,  5.02697766e-01, -1.87914821e+00],
                                  [ 3.24316111e+00, -1.35305546e+00,  5.59683825e-01],
                                  [ 4.49272713e+00,  2.25802114e-01, -8.05304903e-01],
                                  [ 5.74466053e+00,  1.75364395e+00, -1.90801145e+00],
                                  [ 7.09622459e+00, -1.87894658e-10,  6.33205622e-08],
                                  [ 0.00000000e+00,  3.50000000e-01,  1.05300000e+00],
                                  [ 0.00000000e+00, -3.50000000e-01,  1.05300000e+00]]),
              "charge": [ 0, 0, 0, 0, 0, 0, 0.4, -0.4 ],
              "BB bonds": [],
              "SC bonds": np.array([[ 0, 1, 0.315, 5000 ]]),
              "constraints": np.array([[ 1, 2, 0.335],
                                       [ 2, 5, 0.412],
                                       [ 4, 5, 0.293],
                                       [ 1, 4, 0.404 ],
                                       [ 2, 4, 0.470 ]]),
              "BB angles": [],
              "SC angles": np.array([[ 0, 1, 2, 120, 60 ],
                                     [ 0, 1, 4, 130, 60 ]]),
              "BB dihedral": [],
              "SC dihedral": np.array([[ 5, 4, 2, 1, 180, 100]]),
              "BB ref": [ 0 ],
              "BB-SC ref": [ 0, 1 ],
              "BB-SC angle": 100,
              "vsites3": { 6: [ 0.35, 0.00, 1.053 ],
                           7: [ 0.00, 0.35, 1.053 ] },
              "exclusion": [ 6, 7 ],
              "int exclusion": True,
              "vsitesn": [3]
              },
          "martini22p": {
              
              } 
        }
    }

Artificial = {
    # "Z":
    #     { "name": "C15",
    #       "martini22": {
              
    #           },
    #       "marstini22": {
    #           "beads": [ "BB", "SC1", "SC2", "SC3", "SC4", "SC5", "H", "O" ],
    #           "FF name": [ "CPBB", "Nda", "C1", "C1", "C1", "C1", "CPNH", "CPCO" ],
    #           "coords": np.array([[ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00],
    #                               [ 3.49905756e+00,  1.93707215e-08, -2.48815373e-08],
    #                               [ 6.99811512e+00,  0.00000000e+00,  0.00000000e+00],
    #                               [ 1.04971729e+01,  0.00000000e+00,  0.00000000e+00],
    #                               [ 1.39962302e+01,  0.00000000e+00,  0.00000000e+00],
    #                               [ 1.74952878e+01,  0.00000000e+00,  0.00000000e+00],
    #                               [ 0.00000000e+00,  3.50000000e-01,  1.05300000e+00],
    #                               [ 0.00000000e+00, -3.50000000e-01,  1.05300000e+00]]),
    #           "charge": [ 0, 0, 0, 0, 0, 0, 0.4, -0.4 ],
    #           "BB bonds": [],
    #           "SC bonds": np.array([[ 0, 1, 0.47, 1250 ],
    #                                 [ 1, 2, 0.47, 1250 ],
    #                                 [ 2, 3, 0.47, 1250 ],
    #                                 [ 3, 4, 0.47, 1250 ],
    #                                 [ 4, 5, 0.47, 1250 ]]),
    #           "constraints": [],
    #           "BB angles": [],
    #           "SC angles": np.array([[ 0, 1, 2, 180, 25],
    #                                  [ 1, 2, 3, 180, 25],
    #                                  [ 2, 3, 4, 180, 25],
    #                                  [ 3, 4, 5, 180, 25]]),
    #           "BB dihedral": [],
    #           "SC dihedral": [],
    #           "BB ref": [ 0 ],
    #           "BB-SC ref": [ 0, 1 ],
    #           "BB-SC angle": 100,
    #           "vsites3": { 6: [ 0.35, 0.00, 1.053 ],
    #                         7: [ 0.00, 0.35, 1.053 ] },
    #           "exclusion": [ 6, 7 ]
    #           },
    #       "martini22p": {
              
    #           }
    #    }
    }


AA = { **Natural , **Artificial }


"""

IGNORE THIS

def is_integer(n):
    try:
        float(n)
    except ValueError:
        return False
    else:
        return float(n).is_integer()

def TopolInfo( index ):
    res = u.residues[index]
    print('    > name:')
    print(res.resname)
    print('"charmm36": {')
    print('  "beads":','[ "'+'", "'.join(res.atoms.names)+'" ],')
    print('  "FF name":','[ "'+'", "'.join(res.atoms.types)+'" ],')
    print('  "coords":')
    for j, line in enumerate(res.atoms.positions - res.atoms.center_of_geometry()):
        if j==0:
            print("np.array([[",', '.join('{:2.8e}'.format(k) for k in line),'],')
        elif j==len(res.atoms.names)-1:
            print('          [' ,', '.join('{:2.8e}'.format(k) for k in line),']]),')
        else:
            print('          [' ,', '.join('{:2.8e}'.format(k) for k in line),'],')
    print('  "charge":','[ '+', '.join([ str(round(k,4)) for k in res.atoms.charges])+' ],')
    print('  "Bonds":')
    bonds = u.bonds.atomgroup_intersection(res.atoms).indices - res.atoms.indices[0]
    bonds = np.array([line for line in bonds if ( np.array(line) >= 0 ).all() and ( np.array(line) <= (res.atoms.indices[-1]-res.atoms.indices[0]) ).all()])
    for j, line in enumerate(bonds):
        if j ==0:
            print('np.array([[ {} ],'.format(', '.join(['{:2d}'.format(k) for k in line])))
        elif j == len(bonds)-1:
            print('          [ {} ]]),'.format(', '.join(['{:2d}'.format(k) for k in line])))
        else:
            print('          [ {} ],'.format(', '.join(['{:2d}'.format(k) for k in line])))
    print('  "SC bonds": [],')
    print('  "constraints": [],')
    print('  "Angles":')
    angles = u.angles.atomgroup_intersection(res.atoms).indices - res.atoms.indices[0]
    angles = angles[angles[:,1].argsort()]
    angles = np.array([line for line in angles if ( np.array(line) >= 0 ).all() and ( np.array(line) <= (res.atoms.indices[-1]-res.atoms.indices[0]) ).all()])
    for j, line in enumerate(angles):
        if j ==0:
            print('np.array([[ {} ],'.format(', '.join(['{:2d}'.format(k) for k in line])))
        elif j == len(angles)-1:
            print('          [ {} ]]),'.format(', '.join(['{:2d}'.format(k) for k in line])))
        else:
            print('          [ {} ],'.format(', '.join(['{:2d}'.format(k) for k in line])))
    print('  "BB angles": [],')
    print('  "Dihedrals":')
    dihed = u.dihedrals.atomgroup_intersection(res.atoms).indices - res.atoms.indices[0]
    dihed = dihed[dihed[:,1].argsort()]
    dihed = np.array([line for line in dihed if ( np.array(line) >= 0 ).all() and ( np.array(line) <= res.atoms.indices[-1]-res.atoms.indices[0] ).all()])
    for j, line in enumerate(dihed):
        if j ==0:
            print('np.array([[ {} ],'.format(', '.join(['{:2d}'.format(k) for k in line])))
        elif j == len(dihed)-1:
            print('          [ {} ]]),'.format(', '.join(['{:2d}'.format(k) for k in line])))
        else:
            print('          [ {} ],'.format(', '.join(['{:2d}'.format(k) for k in line])))
    print('  "Impropers":')
    d = u.impropers.atomgroup_intersection(res.atoms)
    dihed = [ [ int(d._ags[i][n].index - res.atoms.indices[0]) if  d._ags[i][n].resid == index+1
               else 'P' if d._ags[i][n].resid == index else 'N' if d._ags[i][n].resid == index+2 else "kill"
               for i in range(len(d._ags)) ] for n in range(len(d._ags[0]))]

    dihed = np.array([ [i for i in line] for line in dihed if ( len( set(line) ) > 2 ) ])

    #dihnm = u.impropers.atomgroup_intersection(res.atoms).names
    #dihed = dihed[dihed[:,1].argsort()]
    #dihed = np.array([line for line in dihed if ( np.array(line) >= -1 ).all() and ( np.array(line) <= res.atoms.indices[-1]-res.atoms.indices[0]+1 ).all()])
    for j, line in enumerate(dihed):
        #print( dihnm[j] )
        if j ==0:
            print('np.array([[ {} ],'.format(', '.join(['{:2d}'.format(int(k)) if is_integer(k) else '"{}"'.format(k) for k in line])))
        elif j == len(dihed)-1:
            print('          [ {} ]]),'.format(', '.join(['{:2d}'.format(int(k)) if is_integer(k) else '"{}"'.format(k) for k in line])))
        else:
            print('          [ {} ],'.format(', '.join(['{:2d}'.format(int(k)) if is_integer(k) else '"{}"'.format(k) for k in line])))
    print('  "BB ref":','[ '+', '.join([ str(k) for k in res.atoms.select_atoms('backbone').indices - res.atoms.indices[0] ])+' ],')
    print('  "SC ref":','[ '+', '.join( [ str(i-res.atoms.indices[0]) for i in u.select_atoms('( bonded (group RES and name CA)) and not backbone', RES=res.atoms).indices ] )+' ],')
    print('  "BB-SC ref": [],')
    print('  "BB-SC angle": 0,')
    print('  "vsites3": { },')
    print('  "exclusion": []')
    print('},')
    print(res.resname)

    return


def TopolInfo( index ):
    res = u.residues[index]
    print('    > name:')
    print(res.resname)
    print('"charmm36": {')
    #print('  "beads":','[ "'+'", "'.join(res.atoms.names)+'" ],')
    #print('  "FF name":','[ "'+'", "'.join(res.atoms.types)+'" ],')
    print('  "beads":')
    names = res.atoms.names
    print('[ "'+'", "'.join(names)+'" ],')
    
    print('  "coords":')

    coords = res.atoms.positions - res.atoms.select_atoms('name CA').positions
    
    Angle = ( index%2 * np.pi ) * np.array([0,1,0])
    #coords = R.from_rotvec( Angle ).apply( coords )
    
    #coords[ 1 ] = coords[ 1 ] [::-1]
    #coords[ -1 ] = coords[ -1 ] [::-1]
    coords = AlignWithAxis( coords,
                            coords[res.atoms.select_atoms('name CB').indices[0] - res.atoms.indices[0]], [1,0,0]  )
    #coords[3:-2] = AlignWithAxis( res.atoms.positions[3:-2,:] - res.atoms.select_atoms('name CA').positions,
    #                       res.atoms.select_atoms('name CB').positions ,[1,0,0]  )

    #coords[ 1 ] = coords[ 1 ] [::-1]
    #coords[ -1 ] = coords[ -1 ] [::-1]
    
    for j, line in enumerate(coords):
        if j==0:
            print("np.array([[",', '.join('{:2.8e}'.format(k) for k in line),'],')
        elif j==len(coords)-1:
            print('          [' ,', '.join('{:2.8e}'.format(k) for k in line),']]),')
        else:
            print('          [' ,', '.join('{:2.8e}'.format(k) for k in line),'],')
    
    print('  "BB ref":','[ '+', '.join([ str(k) for k in res.atoms.select_atoms('backbone').indices - res.atoms.indices[0] ])+' ],')
    print('  "SC ref":','[ '+', '.join( [ str(i-res.atoms.indices[0]) for i in u.select_atoms('( bonded (group RES and name CA)) and not backbone', RES=res.atoms).indices ] )+' ],')
    print(res.resname)

    return

def DiffITPs( Mine: str, Other: str):
    for i, field in [(2,'bonds'),(3,'angles'),(4,'dihedrals'),(2,'pairs'),(5,'cmap')]:
        Values_Mine = []
        with open(Mine) as MINE:
            Store = False
            for line in MINE:
                if not line.split():
                    Store = False
                if Store:
                    if not ';' in line.split()[0]:
                        Values_Mine.append( line.split()[:i] )
                if field in line:
                    Store = True
                    
        Values_P2G = []
        with open(Other) as OTHER:
            Store = False
            for line in OTHER:
                if not line.split():
                    Store = False
                if Store:
                    if not ';' in line.split()[0]:
                        Values_P2G.append( line.split()[:i] )
                if field in line:
                    Store = True
                    
        Values = []
        for i in Values_Mine:
            state = 0
            for j in Values_P2G:
                if i == j or i == j[::-1]:
                    state=1
            if not state:
                Values.append( i )
                
        for i in Values_P2G:
            state = 0
            for j in Values_Mine:
                if i == j or i == j[::-1]:
                    state=1
            if not state:
                Values.append( i )
                        
        print( field )
        print( Values )
    return

"""
