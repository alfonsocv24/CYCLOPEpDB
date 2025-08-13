#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 15:55:42 2024

@author: alfonsocabezonvizoso
"""

'''
This code is an adaptation of the CP_Builder.py code developed by Fabián Suárez-Lestón
fabian.suarez.leston@usc.es. This code creates the PDB file of a CP sequence
at the MA(R/S)TINI3 resolution.
'''

import time
import argparse
import subprocess
import MDAnalysis
import numpy as np

from CPDex import AA, FFs

from scipy.spatial.transform import Rotation as R

import warnings
warnings.filterwarnings('ignore')

program_version = '3.0'

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# SCRIPT ARGUMENTS
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# Program description
parser = argparse.ArgumentParser(description =
    '''Cyclopeptide Builder, v{}\n
    Create the PDB and ITP files of a CP or NT given its sequence(s).'''.format( program_version ) )

parser.add_argument( "-s", "--sequence", default = 'A',             #### NEEDED
    help = """Sequence of CPs.""")
parser.add_argument( "-ff", "--forcefield", type = str, default = 'marstini3',
    help = """Forcefield.""")

args = parser.parse_args()


#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# FUNCTIONS
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


def AlignWithAxis( Positions, Vector: np.ndarray, Axis: np.ndarray ):
    '''
    Perform an alignment of a set of points

    Parameters
    ----------
    Positions : np.ndarray
        Coordinates of N points in space. Dimension Nx3.
    Vector : np.ndarray
        Vector to align. Dimension 1x3.
    Axis : np.ndarray
        Reference axis to perform the aligment.

    Returns
    -------
    np.ndarray
        Coordinates of N points in space after the rotation, where Vector
        has become parallel to Axis. Dimension Nx3.

    '''
    
    # Check that Vector and Axis are not parallel
    if np.cross( Vector, Axis ).any():
        # Angle between the vector and the axis
        Angle = np.arccos( np.inner( Vector, Axis ) / ( np.linalg.norm( Vector ) * np.linalg.norm( Axis ) ) )
        
        # Rotation axis, a vector orthogonal to both Vector and Axis of module Angle
        RotAxis = Angle * np.cross( Vector, Axis ) / np.linalg.norm( np.cross( Vector, Axis ) )
        
        # Define the rotation
        Rotation = R.from_rotvec( RotAxis )
        
        return Rotation.apply( Positions )
    
    # If Vector and Axis are antiparallel
    elif np.dot(Vector,Axis) < 0:
        # Rotation axis, a vector orthogonal to both Vector and Axis of module Angle
        RotAxis = np.cross( np.cross( Vector, Vector + np.array([0,0,1])), Vector )
        RotAxis *= np.pi/np.linalg.norm( RotAxis )
        
        # Define the rotation
        Rotation = R.from_rotvec( RotAxis )
        
        return Rotation.apply( Positions )
        
    # If Vector and Axis are already parallel
    else:
        return Positions



class Cyclopeptide:
    '''
    Generation of a cyclopeptide
    
    Attributes
    ----------
    name : str
        A label to identify the CP
    sequence : str
        The sequence of amino acids of the CP, starting in a D-amino acid.
    ff : str
        A forcefield for the structure of the CP.
    length: int
        The number of residues in the CP. Is the lenght of the sequence.
        
    Methods
    -------
    CheckSequence
        Check if the sequence of amino acids is valid.
    CPTemplate
        Generate a template for the distribution of amino acids.
    GenerateGeometry
        Build a initial geometry of the CP.
    GenerateTopology
        Build the topology of the CP for a simulation with using the given FF.
    OptimizeStructure
        Perform a energy minimization to obtain a more stable structure.
        
    '''
    
    def __init__( self, sequence, forcefield ):
        '''
        Defines the atributes of the CP

        Parameters
        ----------
        sequence : str
            The sequence of amino acids of the CP, starting in a D-amino acid.
        forcefield : str
            A forcefield for the structure of the CP.
            
        '''
        
        self.sequence = sequence
        self.ff = forcefield
        self.length = len( sequence )
    
    
    def __str__( self ):
        return f'Cyclopeptide of length {self.length} in the {self.ff} FF. Sequence: {self.sequence}'
    
    
    def CheckSequence( self ):
        '''
        Check if the sequence of amino acids is valid.

        Returns
        -------
        bool
            Validity of the sequence.

        '''
        
        for AminoAcid in self.sequence:
            if AminoAcid not in AA.keys():
                return False
        return True
    
    
    def CPTemplate( self ) -> np.ndarray:
        '''
        Generate a template for the distribution of amino acids.

        Returns
        -------
        BB : np.ndarray
            Coordinates of the reference elements of the backbone.

        '''
        
        # Key of the forcefield ("AA" if not a CG FF)
        FF = self.ff if ( self.ff in FFs.keys() and FFs[self.ff]["CG"]=='yes' ) else "AA"
        
        # Number of elements in the backbone for each amino acid
        AASpace = np.array( [ len( AA[AminoAcid][FF]["BB ref"] ) for AminoAcid in self.sequence ] )
        
        # Number of elements in the backbone
        n = np.sum( AASpace )
        
        # Generate n evenly spaced points in a circle of radius l/2sin(π/n)
        BBPos = np.array([ [FFs[FF]["distance"]/(2*np.sin(np.pi/n))*np.cos(i), FFs[FF]["distance"]/(2*np.sin(np.pi/n))*np.sin(i),0]
                          for i in np.linspace( 0, 2*np.pi, n, endpoint=False ) ])
        
        # Average the position of the positions for larger residues
        BB=np.array( [ np.mean( BBPos[int(np.sum(AASpace[:i])):BBl+int(np.sum(AASpace[:i])),:],axis=0 )
                       for i, BBl in enumerate( AASpace ) ] )
        
        return BB
    
    
    def GenerateGeometry( self ):
        '''
        Build a initial geometry of the CP.

        Returns
        -------
        None.
            Writes a PDB file with the structure.

        '''
        
        # A initial template for the distribution of amino acids
        BB = self.CPTemplate()
        
        # Key of the forcefield ("AA" if not a CG FF)
        FF = self.ff if ( self.ff in FFs.keys() and FFs[self.ff]["CG"]=='yes' ) else "AA"

        # Initialize an array for storing the coordinates of the CP
        CPCoords = np.array([0,0,0])
        
        # Initialize a phase for the alternance between D and L-amino acids
        #   It must be modified if later on other patterns due to the adition
        #   of non-natural α-amino acids is incorporated.
        Phase = 1
        
        # Iterate over the elements of the CP
        for i, AminoAcid in enumerate( self.sequence ):
            # The angle between the initial position (accoding to the template) and the X-axis
            Angle = np.arccos( np.inner( BB[i,:], np.array([1,0,0]) ) / np.linalg.norm( BB[i,:] ) ) 
            
            # The angle-axis vector for a rotation along the Z-axis
            Angle *= np.array([0,0,1])
            
            # Ensure a left-hand rotation by changing the sign of the angle if required
            Angle = -Angle if np.dot( np.cross( np.array([1,0,0]), BB[i,:] ), np.array([0,0,1]) ) < 0 else Angle
            
            # Incorporate the coordinates of the residue, properly aligned and with the right chirality
            CPCoords = np.vstack( [ CPCoords, *( R.from_rotvec( Angle ).apply( 
                np.hstack( [AA[AminoAcid][FF]["coords"][:,:2], Phase * AA[AminoAcid][FF]["coords"][:,2:] ] )
                ) + BB[i] ) ] )
            
            # Change the value of the Phase to alternate the chirality of the next residue
            #   It must be modified if later on other patterns due to the adition
            #   of non-natural α-amino acids is incorporated.
            Phase *= -1
            
        # Remove the first (dummy) line of the coordinates
        CPCoords = CPCoords[1:]

        # An array with the index, name of the residue and name of the bead 
        CPInfo = np.array([ [i, AA[AminoAcid]["name"], bead ] for i, AminoAcid in enumerate(self.sequence) 
                           for bead in AA[AminoAcid][FF]["beads"] ])
        
        # Define a MDAnalysis.Universe object with the information of the CP
        CPStructure = MDAnalysis.Universe.empty( len(CPInfo), 
                                                n_residues=self.length, 
                                                atom_resindex=CPInfo[:,0],
                                                residue_segindex=[0]*self.length,
                                                trajectory=True )
        
        # Add the atributes and coodinates
        CPStructure.add_TopologyAttr( 'resname', [ AA[i]["name"] for i in self.sequence ] )
        CPStructure.add_TopologyAttr( 'name', CPInfo[:,2] )
        CPStructure.add_TopologyAttr( 'resid', [ i+1 for i in range(len(self.sequence)) ] )
        CPStructure.atoms.positions = CPCoords
        
        # Write into an extenral file
        CPStructure.atoms.write( f'{self.sequence}.pdb' )
        
        return
        
    
    def GenerateTopology( self ):
        '''
        Build the topology of the CP for a simulation with using the given FF.

        Returns
        -------
        None.
            Writes a ITP file with the structure.

        '''
        global AminoIndex
        global BondIndices, AngleIndices, DihedralIndices
        Backbone = []; AminoIndex = []
        
        # Key of the forcefield ("AA" if not a CG FF)
        FF = self.ff if ( self.ff in FFs.keys() and FFs[self.ff]["CG"]=='yes' ) else "AA"
        
        # Generate the topology of a CG peptide
        if ( self.ff in FFs.keys() and FFs[self.ff]["CG"]=='yes' ):
            with open( '{}.itp'.format( self.sequence ), 'w+' ) as ITP:
    ### Header
                ITP.write("; Topology file for the ciclopeptide:\n")
                ITP.write("; {}\n".format( self.sequence ))
                ITP.write("; Sequence of the cyclopeptide:\n")
                ITP.write("; {}\n".format( self.sequence ))
                ITP.write("; Forcefield:\n")
                ITP.write("; {}\n".format( self.ff ))
                ITP.write("; File generated on {} using the CPBuilder v{}\n".format( time.strftime("%d/%m/%Y %H:%M:%S"), program_version ))
                ITP.write("; If you use this file in your research, please, include the corresponding citation\n\n")
                
    ### [ moleculetype ]
                ITP.write("[ moleculetype ]\n")
                ITP.write("; Name          Number\n")
                ITP.write("{}    1\n\n".format(self.sequence))
                
    ### [ atoms ]
                ITP.write("[ atoms ]\n")
                
                BeadCount = 1
                # Iterate over the residues of the CP
                for i, AminoAcid in enumerate( self.sequence ):
                    Amino = AA[AminoAcid][self.ff]
                    
                    # Iterate over the beads of the AA
                    for j in range(len(Amino["beads"])):
                        # Store the beads of the backbone in other list
                        if j in Amino["BB ref"]:
                            Backbone.append( BeadCount )
                        
                        # Store the index of the first element of a new residue
                        if j == 0:
                            AminoIndex.append( BeadCount )
                            
                        # Write a line in the [ atoms ] section
                        ITP.write( ('{:>5} {:<6} {:>3}' + ' {:>6}'*3 + ' {:>7.4f} {} \n').format(
                            BeadCount,                                         # Index of the bead, 1-leading
                            Amino["FF name"][j],                               # Name of the bead in the FF
                            i+1,                                               # Index of the residue, 1-leading
                            AA[AminoAcid]["name"],                             # Name of the residue
                            Amino["beads"][j],                                 # Name of the bead
                            i+1,                                               # Index of the residue, 1-leading
                            Amino["charge"][j],                                # Charge of the bead
                            '{:>7.4f}'.format(Amino["mass"][j]) if Amino["mass"][j] != None else "" # Mass of the bead, if definied
                        ) )
                        
                        # Update the bead counter
                        BeadCount += 1
                        
                ITP.write("\n")
                
    ### [ bonds ]
                ITP.write("[ bonds ]\n")
                
    # (CG) Backbone bonds
                ITP.write("; Bonds between elements of the backbone\n")
                
                for i in range(len(Backbone)):
                    ITP.write("{:>5} {:>6} {:>7} {:>10.5f} {:>6} ; {:<14}\n".format(
                            Backbone[i],                                       # Index of the i-th element of the backbone
                            Backbone[(i+1)%len(self.sequence)],                # Index of the (i+1)-th element of the backbone
                            1,                                                 # Bond function
                            0.38,                                              # Bond length
                            2275 if ( self.ff == 'marstini3' or self.ff == 'martini3' ) else 6275, # Force constant, maybe better define it in CPDex
                        '{0}(BB) - {1}(BB)'.format(
                            AA[self.sequence[i]]["name"],
                            AA[self.sequence[(i+1)%len(self.sequence)]]["name"]
                            )                                                  # A comment describing the residues being bonded by this entry
                        ))

    # (CG) Backbone-sidechain & sidechain-sidechain bonds
                ITP.write("; Bonds between elements of the same residue\n")
                
                for i, AminoAcid in enumerate( self.sequence ):
                    Bonds = AA[AminoAcid][self.ff]["SC bonds"]
                    
                    for bond in Bonds:
                        ITP.write("{:>5.0f} {:>6.0f} {:>7} {:>10.5f} {:>6.0f} ; {:<14}\n".format(
                                *bond[:2] + AminoIndex[i],                     # Index of the i-th and (i+1)-th elements of the bond
                                1,                                             # Bond function
                                *bond[2:],                                     # Bond lenght and force constant
                                f'{AA[AminoAcid]["name"]} internal bond'       # A comment describing the elements being bonded by this entry
                            ))
                        
                
                
    # (CG) Sidechain-sidechain bonds (only applied when constraints appear)
                ITP.write("; Bonds associated to constraints\n")
                ITP.write("#ifdef FLEXIBLE\n")
                for i, AminoAcid in enumerate( self.sequence ):
                    Constraints = AA[AminoAcid][self.ff]["constraints"]
                    
                    for constraint in Constraints:
                        ITP.write("{:>5.0f} {:>6.0f} {:>7} {:>10.5f} {:>10d} ; {:<14}\n".format(
                                *constraint[:2] + AminoIndex[i],               # Index of the i-th and (i+1)-th elements of the constraint/bond
                                1,                                             # Bond function
                                constraint[2],                                 # Bond length
                                1000000,                                       # Force constant
                                AA[AminoAcid]["name"]                          # A comment with the name of the amino acid
                            ))
                
                ITP.write("#endif\n\n")
                
    ### [ constraints ]
                ITP.write("[ constraints ]\n")
                ITP.write("#ifndef FLEXIBLE\n")
                for i, AminoAcid in enumerate( self.sequence ):
                    Constraints = AA[AminoAcid][self.ff]["constraints"]
                    
                    for constraint in Constraints:
                        ITP.write("{:>5.0f} {:>6.0f} {:>7} {:>10.5f} ; {:<14}\n".format(
                                *constraint[:2] + AminoIndex[i],               # Index of the i-th and (i+1)-th elements of the constraint
                                1,                                             # Constraint function
                                *constraint[2:],                               # Constraint length
                                AA[AminoAcid]["name"]                          # A comment with the name of the amino acid
                            ))
                        
                ITP.write("#endif\n\n")
                
    ### [ virtual_sites3 ]
                ITP.write("[ virtual_sites3 ]\n")
                
                Sign = -1
                for i, AminoAcid in enumerate(self.sequence):
                    VSites = AA[AminoAcid][self.ff]["vsites3"]
                    
                    for site in VSites:
                        ITP.write( '{:>5} {:>6} {:>6} {:>6} {:>6} {:>7.4f} {:>7.4f} {:>7.4f} ; {:<14}\n'.format(
                            site + AminoIndex[i],                              # Index of the virtual site
                            Backbone[i],                                       # Index of the BB bead it is linked to
                            Backbone[i-1],                                     # Index of the next BB bead
                            Backbone[(i+1)%len(self.sequence)],                # Index of the previous BB bead
                            4,                                                 # Type of virtual site
                            *VSites[site][:2],                                 # X and Y coordinates of the virtual site
                            Sign * VSites[site][2],                            # Z coordinate of the virtual site
                            "{0}({1})".format(AA[AminoAcid]["name"],AA[AminoAcid][self.ff]["beads"][site]) # A comment with information about the virtual site
                            ))
                        
                    # Change the sign of the Z coordinate of the virtual site to stablish D and L-amino acid alternance
                    Sign = -Sign
                    
                ITP.write("\n")
                
    ### [ virtual_sitesn ]
                ITP.write("[ virtual_sitesn ]\n")
                
                for i, AminoAcid in enumerate(self.sequence):
                    for site in AA[AminoAcid][self.ff]["vsitesn"]:
                        ITP.write( '{:>5} {:>6} {:>6} {:>6} {:>6} {:>6}\n'.format(
                            site + AminoIndex[i],                              # Index of the virtual site
                            2,                                                 # Type of virtual site
                            site + AminoIndex[i] + 2,                          # Index of the bead two positions forward
                            site + AminoIndex[i] + 1,                          # Index of the bead one positions forward
                            site + AminoIndex[i] - 1,                          # Index of the bead one positions backwards
                            site + AminoIndex[i] - 2                           # Index of the bead two positions backwards
                            ))
                        
                ITP.write("\n")

                
    ### [ exclusions ]
                ITP.write("[ exclusions ]\n")
                
                for i, AminoAcid in enumerate(self.sequence):
                    # Internally defined exclusions
                    if AA[AminoAcid][self.ff]["int exclusion"]:  ### CHECK!
                        # Generate exclusions between elements in the side chain
                        for j in range(len(AA[AminoAcid][self.ff]["beads"])-(3 if "marstini" in self.ff else 1)):
                            ITP.write( " ".join( list( map( str, AminoIndex[i] + np.arange(j,len(AA[AminoAcid][self.ff]["beads"])-(2 if "marstini" in self.ff else 0),1) ) ) ) + "\n")
                    
                    # Exclusions corresponding to the virtual sites
                    Exclusions = VSites = AA[AminoAcid][self.ff]["exclusion"]
                    
                    for excl in Exclusions:
                        ITP.write( '{:>5}  {}\n'.format(
                            excl + AminoIndex[i],                              # The index of the bead in the exclusion
                            '  '.join([str(j+1) for j in range(BeadCount-1) if j+1 != (excl+AminoIndex[i])])  # The index of the rest of the beads in the CP
                        ))
                        
                ITP.write("\n")
                
    ### [ angle ]
                ITP.write("[ angles ]\n")
                
    # (CG) Backbone angles
                ITP.write("; Angles between elements of the backbone\n")
                
                for i in range(len(Backbone)):
                    ITP.write( '{:>5} {:>6} {:>6} {:>7} {:>7} {:>6} ; {:<14}\n'.format( 
                        Backbone[i],                                           # Index of the i-th bead of the backbone
                        Backbone[(i+1)%len(self.sequence)],                    # Index of the (i+1)-th bead of the backbone
                        Backbone[(i+2)%len(self.sequence)],                    # Index of the (i+2)-th bead of the backbone
                        1,                                                     # Angle function
                        (len(self.sequence) - 2)/len(self.sequence)*180,       # Value of the angle
                        227 if ( self.ff == 'marstini3' or self.ff == 'martini3' ) else 627, # Force constant
                        '{0}(BB) - {1}(BB) - {2}(BB)'.format(
                            AA[self.sequence[i]]["name"],
                            AA[self.sequence[(i+1)%len(self.sequence)]]["name"],
                            AA[self.sequence[(i+2)%len(self.sequence)]]["name"]
                            )                                                  # A comment describing the angle
                       ))

    # (CG) Backbone - sidechain angles
                ITP.write("; Angles between elements of the backbone and side chains\n")

                for i, AminoAcid in enumerate( self.sequence ):
                    Angle = AA[AminoAcid][self.ff]["BB-SC ref"]
                    
                    if Angle:
                        ITP.write( '{:>5} {:>6} {:>6} {:>7} {:>7} {:>6} ; {:<14}\n'.format( 
                            Angle[1] + AminoIndex[i],                          # Index of the first element of the side chain
                            Angle[0] + AminoIndex[i],                          # Index of the i-th bead of the backbone
                            Backbone[(i+1)%len(self.sequence)],                # Index of the (i+1)-th bead of the backbone
                            2,                                                 # Angle function
                            (360-((len(self.sequence)-2)*180/len(self.sequence)))/2, # Value of the angle
                            25,                                                # Force constant
                            '{0}(BB) - {0}(BB) - {1}(SC)'.format(
                                AA[self.sequence[i-1]]["name"], 
                                AA[self.sequence[i]]["name"]
                            )                                                  # A comment describing the angle
                        ))

    # (CG) Sidechain angles
                ITP.write("; Angles within the side chains\n")
                
                for i, AminoAcid in enumerate( self.sequence ):
                    Angles = AA[AminoAcid][self.ff]["SC angles"]
                    
                    for angle in Angles:
                        ITP.write( '{:>5} {:>6} {:>6} {:>7} {:>7} {:>6} ; {:<14}\n'.format( 
                            *angle[:3] + AminoIndex[i],                        # Indices of the elements defining the angle
                            2,                                                 # Angle function
                            *angle[3:],                                        # Value of the angle and force constant
                            AA[self.sequence[i]]["name"]                       # A comment describing the angle
                        ))
                        
                ITP.write("\n")
                
    ### [ dihedrals ]
                ITP.write("[ dihedrals ]\n")
                
    # (CG) Backbone dihedrals
                ITP.write("; Dihedrals between elements of the backbone\n")
                for i in range(len(Backbone)):
                    ITP.write( '{:>5} {:>6} {:>6} {:>6} {:>7} {:>7} {:>6} ; {:<}\n'.format( 
                        Backbone[i],                                           # Index of the i-th bead of the dihedral
                        Backbone[(i+1)%len(self.sequence)],                    # Index of the (i+1)-th bead of the dihedral
                        Backbone[(i+2)%len(self.sequence)],                    # Index of the (i+2)-th bead of the dihedral
                        Backbone[(i+3)%len(self.sequence)],                    # Index of the (i+3)-th bead of the dihedral
                        2,                                                     # Dihedral function
                        0,                                                     # Value of the dihedral angle
                        100 if ( self.ff == 'marstini3' or self.ff == 'martini3' ) else 418, # Force constant       
                        '{0}(BB) - {1}(BB) - {2}(BB) - {3}(BB)'.format(
                            AA[self.sequence[i]]["name"],
                            AA[self.sequence[(i+1)%len(self.sequence)]]["name"],
                            AA[self.sequence[(i+2)%len(self.sequence)]]["name"],
                            AA[self.sequence[(i+3)%len(self.sequence)]]["name"]
                            )                                                  # A comment describing the dihedral
                        ))
                    
    # (CG) Backbone - sidechain dihedrals
                if self.ff == 'marstini3' or self.ff == 'martini3':  ## Fonso
                    pass
                
                else:
                    ITP.write("; Dihedrals involving both elements of the backbone and side chains\n")
                    for i, AminoAcid in enumerate( self.sequence ):
                        SC = AA[AminoAcid][self.ff]["BB-SC ref"]
                        
                        if SC:
                            ITP.write( '{:>5} {:>6} {:>6} {:>6} {:>7} {:>7} {:>6} ; {:<14}\n'.format( 
                                Backbone[i-1],                                 # Index of the (i-1)-th bead of the backbone
                                SC[0] + AminoIndex[i],                         # Index of the i-th bead of the backbone
                                Backbone[(i+1)%len(self.sequence)],            # Index of the (i+1)-th bead of the backbone
                                SC[1] + AminoIndex[i],                         # Index of the first element of the side chain
                                2,                                             # Dihedral function
                                180,                                           # Value of the dihedral angle
                                418,                                           # Force constant
                                '{0}(BB) - {1}(BB) - {2}(BB) - {1}(SC)'.format(
                                    AA[self.sequence[i-1]]["name"], 
                                    AA[self.sequence[i]]["name"], 
                                    AA[self.sequence[(i+1)%len(self.sequence)]]["name"]
                                )                                              # A comment describing the dihedral
                            ))
                            
    # (CG) Sidechain dihedrals
                ITP.write("; Dihedrals within the side chains\n")
                for i, AminoAcid in enumerate( self.sequence ):
                    SC = AA[AminoAcid][self.ff]["SC dihedral"]
                    
                    for dihedral in SC:
                        ITP.write( '{:>5} {:>6} {:>6} {:>6} {:>7} {:>7} {:>6} ; {:<14}\n'.format(
                            *dihedral[:4] + AminoIndex[i],                     # Indices of the elements in the sidechain defining a dihedral
                            2,                                                 # Dihedral function
                            *dihedral[4:],                                     # Value of the dihedral angle and force constant
                            AA[AminoAcid]["name"]                              # A comment identifying the dihedral
                        ))
                        
                ITP.write("\n")
                
    # ### [ cmap ]  Why in hell is this section here!?
    #             ITP.write("[ cmap ]\n")
    #             if FFs[self.ff]['CG'] == 'no':
    #                 for i, amino in enumerate(self.sequence):
    #                     ITP.write( '{:>5} {:>6} {:>6} {:>6} {:>6} {:>7} ; {:<14}\n'.format(
    #                         AA[self.sequence[i-1]][self.ff]["BB ref"][-2] + AminoIndex[i-1],
    #                         *[ j + AminoIndex[i] for j in AA[self.sequence[i]][self.ff]["BB ref"][:-1] ],
    #                         AA[self.sequence[(i+1)%len(self.sequence)]][self.ff]["BB ref"][0] + AminoIndex[(i+1)%len(self.sequence)],
    #                         1,
    #                         ''
    #                         ))
    #             ITP.write("\n")
    
        # Generate the topology of an AA/UA peptide
        else:
            # Iterate over the sequence of AAs
            BeadCount = 1
            for i, AminoAcid in enumerate( self.sequence ):
                Amino = AA[AminoAcid][FF]
                
                # Iterate over the beads of the AA
                for j in range(len(Amino["beads"])):
                    # Store the index of the atoms in the backbone 
                    if j in Amino["BB ref"]:
                        Backbone.append( BeadCount )
                        
                    # Store the index of the first element of the residue
                    if j == 0:
                        AminoIndex.append( BeadCount )
                        
                    # Update the counter of beads
                    BeadCount += 1
        
            # Execute pdb2gmx in order to generate a valid topology of the CP
            subprocess.run(
                'gmx pdb2gmx -f {} -o {} -p {} -ff {} -water none -ignh'.format(
                    '{}.pdb'.format( self.sequence ),              # Structure file of the CP
                    '{}.pdb'.format( self.sequence ),              # New structure, updates the old file
                    '{}.top'.format( self.sequence ),              # Topology of the CP in the FF
                    self.ff,                                                   # Name of the FF
                    ),
                shell = True,
                stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
            
            # Load the structure of the CP
            u = MDAnalysis.Universe( f'{self.sequence}.pdb' )
            
            # Place N-H and C-O bonds perpendicular to the plane of the CP
            Phase = 1
            for i, AminoAcid in enumerate( self.sequence ):
                # The coordinates of N and H 
                N  = u.select_atoms('resid {} and name N'.format(i+1)).positions
                Hn = u.select_atoms('resid {} and (name HN or name H)'.format(i+1)).positions
                
                # Distance between N and H
                DNH = np.linalg.norm( Hn - N )
                
                # Update the position of the Hn to point in the Z direction
                u.select_atoms('resid {} and name HN'.format(i+1)).positions = u.select_atoms('resid {} and name N'.format(i+1)).positions + Phase * np.array([0,0,DNH])

                # The coordinates of C and O
                C  = u.select_atoms('resid {} and name C'.format(i+1)).positions
                O = u.select_atoms('resid {} and name O'.format(i+1)).positions
                
                # Distance between C and O
                DCO = np.linalg.norm( O - C )
                
                # Update the position of the O to point in the Z direction
                u.select_atoms('resid {} and name O'.format(i+1)).positions = u.select_atoms('resid {} and name C'.format(i+1)).positions + Phase * np.array([0,0,DCO])
            
                # Change the phase to allow the D and L-amino acid alternance
                Phase *= -1
                
            # Rewrite the structure file
            u.atoms.write( f'{self.sequence}.pdb' )
            
            Lines = []; Read = False
            with open( f'{self.sequence}.top', 'r+' ) as TOP:
                for line in TOP:
                    # Once the field [ atoms ] appears, start to store the data
                    if "[ atoms ]" in line:
                        Read = True
                    # Once the field [ system ] appears, stop reading data
                    if "[ system ]" in line:
                        break
                    # Read the content of the topology file between the []
                    if Read:
                        Lines.append( line )
            
            # Create a new topology file with ITP extension
            with open( f'{self.sequence}.itp', 'w+' ) as ITP:
                # Add a [ moleculetype ] section in the beginning
                ITP.write("[ moleculetype ]\n")
                ITP.write("; Name          Number\n")
                ITP.write(f"{self.sequence}    1\n\n")
                
                # Store the data recovered from the TOP file
                for line in Lines:
                    ITP.write( line )
            
            # Generate a list with the index each residue begins in
            AminoIndex = [ res.atoms.indices[0]+1 for res in MDAnalysis.Universe('{}.pdb'.format( self.sequence ) ).residues ]
            
            # Remove the TOP file
            subprocess.run( f'rm {self.sequence}.top',
                shell = True )
        
        return
            
        
    def OptimizeStructure( self ):
        '''
        Perform a energy minimization to obtain a more stable structure.

        Returns
        -------
        None.

        '''
        
        # Key of the forcefield ("AA" if not a CG FF)
        FF = self.ff if ( self.ff in FFs.keys() and FFs[self.ff]["CG"]=='yes' ) else "AA"

        # Create a file to add restraints to the structure
        with open( 'restraints.itp', "w+" ) as TOP:
            TOP.write( '[ dihedral_restraints ]\n')
            if FF == 'AA':
                for i, amino in enumerate(self.sequence):
                    # Honestly, I don't know what this criterium means; but it must be working properly. Blame my past self.
                    if len( set( [ AA[self.sequence[i-1]][FF]["BB ref"][-1] + AminoIndex[i-1],
                                   AA[self.sequence[i-1]][FF]["BB ref"][2] + AminoIndex[i-1],
                                   AA[self.sequence[i]][FF]["BB ref"][0] + AminoIndex[i],
                                   AA[self.sequence[i]][FF]["BB ref"][0] + AminoIndex[i] + 1 ] ) ) == 4:
                    # O(i-1)-C(i-1)-N(i)-Hn(i)
                        TOP.write( '{:>5} {:>6} {:>6} {:>6} {:>6} {:>6} {:>6} {:>7} ; {:<14}\n'.format( 
                            AA[self.sequence[i-1]][FF]["BB ref"][-1] + AminoIndex[i-1],  # Index of atom O(i-1)
                            AA[self.sequence[i-1]][FF]["BB ref"][2] + AminoIndex[i-1],   # Index of atom C(i-1)
                            AA[self.sequence[i]][FF]["BB ref"][0] + AminoIndex[i],       # Index of atom N(i)
                            AA[self.sequence[i]][FF]["BB ref"][0] + AminoIndex[i] + 1,   # Index of atom Hn(i). If Hn is not direcly after N it won't work!
                            1,                                                           # Dihedral function
                            180,                                                         # Dihedral angle
                            0,                                                           # Angular threshold for activation of the dihedral restraint  
                            150,                                                         # Force constant
                            ''
                            ))
                        
                    # Hn(i)-N(i)-C(i)-O(i)
                        TOP.write( '{:>5} {:>6} {:>6} {:>6} {:>6} {:>6} {:>6} {:>7} ; {:<14}\n'.format( 
                            AA[self.sequence[i]][FF]["BB ref"][0] + AminoIndex[i] + 1,   # Index of atom Hn(i). If Hn is not direcly after N it won't work!
                            AA[self.sequence[i]][FF]["BB ref"][0] + AminoIndex[i],       # Index of atom N(i)
                            AA[self.sequence[i]][FF]["BB ref"][2] + AminoIndex[i],       # Index of atom C(i-1)
                            AA[self.sequence[i]][FF]["BB ref"][-1] + AminoIndex[i],      # Index of atom O(i-1)
                            1,                                                           # Dihedral function
                            0,                                                           # Dihedral angle
                            0,                                                           # Angular threshold for activation of the dihedral restraint 
                            150,                                                         # Force constant
                            ''
                            ))
                        
                    # Hn(i-2)-N(i-2)-C(i)-O(i)
                        TOP.write( '{:>5} {:>6} {:>6} {:>6} {:>6} {:>6} {:>6} {:>7} ; {:<14}\n'.format( 
                            AA[self.sequence[i]][FF]["BB ref"][0] + AminoIndex[i-2] + 1, # Index of atom Hn(i-2). If Hn is not direcly after N it won't work!
                            AA[self.sequence[i]][FF]["BB ref"][0] + AminoIndex[i-2],     # Index of atom N(i-2)
                            AA[self.sequence[i]][FF]["BB ref"][2] + AminoIndex[i],       # Index of atom C(i)
                            AA[self.sequence[i]][FF]["BB ref"][-1] + AminoIndex[i],      # Index of atom O(i)
                            1,                                                           # Dihedral function
                            0,                                                           # Dihedral angle
                            0,                                                           # Angular threshold for activation of the dihedral restraint 
                            150,                                                         # Force constant
                            ''
                            ))

        TOP.close()
        
        # Generate a new TOP file for the CP
        with open( 'topol.top', "w+" ) as TOP:
            TOP.write( '#include "{}"\n'.format( 
                './{}.itp'.format(self.ff) if FFs[FF]['where']=='itp' 
                else '{}.ff/forcefield.itp'.format(self.ff) if FFs[FF]['where']=='defined'
                else 'NULL' ) )
            TOP.write( '#include "./{}.itp"\n'.format( self.sequence ) )
            TOP.write( '#include "./restraints.itp"\n')
            TOP.write( '\n' )
            TOP.write( '[ system ]\n')
            TOP.write( 'Cyclation of {}\n'.format( self.sequence ) )
            TOP.write( '\n' )
            TOP.write( '[  molecules ]\n' )
            TOP.write( ';name              number\n' )
            TOP.write( '{}       1\n'.format( self.sequence ) )
            TOP.write( '\n' )
        TOP.close()
        
        # # If the FF is in an external file, copy it into the directory
        # if FFs[FF]['where']=='itp' :
        #     subprocess.run(
        #         'cp ./{}.itp .'.format(
        #             self.ff                                                  
        #             ),
        #         shell = True )

        
        # Generation of the box for the minimization in vacuo
        subprocess.run(
            'gmx editconf -f {} -o {} -box {}'.format(
                f'{self.sequence}.pdb',                              # Peptide file
                'peptide.gro',                                                 # Peptide-in-box file
                ' '.join([ str(20*0.38/np.sin(np.pi/len(self.sequence))) ]*3)  # Dimensions of the box
                ),
            shell = True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)

        # Generation of the input file for the minimization in vacuo
        subprocess.run(
            'gmx grompp -p {} -c {} -f {} -o {} -maxwarn 1'.format(
                'topol.top',                                                   # Topology file
                'peptide.gro' ,                                                # Peptide file
                'steep.mdp',                             # MDP file for the minimization
                'minimization.tpr'                                             # Executable for the minimization
                ),
            shell = True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True )
        
        # Minimization of the energy of the system
        subprocess.run(
            'gmx mdrun -deffnm minimization ', 
            shell = True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
        
        # Buiding of the index file
        subprocess.run( ' echo "q" | gmx make_ndx -f minimization -o index.ndx',
            shell = True , stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
        
        # Process the coordinates and center in the box 
        u = MDAnalysis.Universe( 'minimization.gro' )
        u.atoms.positions -= u.select_atoms( FFs[FF]["backbone"] ).center_of_geometry()
        
        if not ( self.ff in FFs.keys() and FFs[self.ff]["CG"]=='yes' ):
            # Place N-H and C-O bonds perpendicular to the plane of the CP
            Phase = 1
            for i, AminoAcid in enumerate( self.sequence ):
                # The coordinates of N and H 
                N  = u.select_atoms('resid {} and name N'.format(i+1)).positions
                Hn = u.select_atoms('resid {} and (name HN or name H)'.format(i+1)).positions
                
                # Distance between N and H
                DNH = np.linalg.norm( Hn - N )
                
                # Update the position of the Hn to point in the Z direction
                u.select_atoms('resid {} and name HN'.format(i+1)).positions = u.select_atoms('resid {} and name N'.format(i+1)).positions + Phase * np.array([0,0,DNH])

                # The coordinates of C and O
                C  = u.select_atoms('resid {} and name C'.format(i+1)).positions
                O = u.select_atoms('resid {} and name O'.format(i+1)).positions
                
                # Distance between C and O
                DCO = np.linalg.norm( O - C )
                
                # Update the position of the O to point in the Z direction
                u.select_atoms('resid {} and name O'.format(i+1)).positions = u.select_atoms('resid {} and name C'.format(i+1)).positions + Phase * np.array([0,0,DCO])
            
                # Change the phase to allow the D and L-amino acid alternance
                Phase *= -1
           
        # Store the modified coordinates
        u.atoms.write( f'{self.sequence}.pdb' )
        
        # Remove some files
        subprocess.run(
            'rm minimization.edr minimization.trr minimization.gro minimization.tpr peptide.gro mdout.mdp minimization.log index.ndx restraints.itp topol.top',
            shell = True )
        
        return
    
    

      
        
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# EXECUTION
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

if __name__ == "__main__":

    # Set the timer
    start_script_time=time.time()
    print(args.sequence)
    
    # Check that the sequence consists in an even number of amino acids
    if len( args.sequence ) % 2 != 0:
        raise ValueError('Only sequences with an even number of amino acids are allowed!')

    # Generate a object for the nanotube
    CP = Cyclopeptide(args.sequence, args.forcefield)
    
    # Assemble the CPs
    CP.GenerateGeometry()
    
    # Generate the nanotube
    CP.GenerateTopology()
    
    CP.OptimizeStructure()
    
    print( f"--- TOTAL TIME: {time.time()-start_script_time:.3f} seconds ---" )
        
