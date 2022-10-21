# -*- coding: utf-8 -*-
"""
Created on Tue Oct 18 11:11:54 2022

@author: 983045
"""
from functions import *
import smiles_syntax

mol = Chem.MolFromSmiles("c1(C=O)cc(OC)c(O)cc1")
#mol = Chem.MolFromSmiles("CC(=O)NC1=C(C=C(C=C1)O)O")
new_mol = cut_bonds_to_core(mol)
print(new_mol)
'''
Above is complete first part
see functions.py
'''
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem.Draw import MolDrawing, DrawingOptions

DrawingOptions.bondLineWidth=1.8
DrawingOptions.atomLabelFontSize=14
DrawingOptions.includeAtomNumbers=False

def get_cores_and_rs(mol):
    main = []
    cores = []
    r_groups = []
    new_mol = cut_bonds_to_core(mol)
    new_mol = new_mol.replace('*', '')
    new_mol = new_mol.replace('()', '')
    new_mol = new_mol.split('.')
    
    for thing in new_mol:
        b = Chem.MolFromSmiles(thing)
        c = has_rings(b)
        if len(c) == 0:
            r_groups.append(thing)
        if len(c) != 0:
            cores.append(thing)
    main.append(cores)
    main.append(r_groups)
    return(main)
new_mol = cut_bonds_to_core(mol)
f = get_cores_and_rs(mol)
'''   
This function takes an rdkit mol object and returns a 2d array with a core 
and a list of r groups in SMILES format
'''
atom_set = get_atom_set(mol) # Gets all atoms in the molceule
atom_ring_set = get_ring_atoms(mol)

atom_set.symmetric_difference_update(atom_ring_set) # creates a set without the ring atoms

core = Chem.MolFromSmiles(f[0][0])
new_mol = cut_bonds_to_core(mol)

boiled = smiles_syntax.convert_wildcards_to_closures(new_mol)
#fried = smiles_syntax.convert_wildcards_to_closures(new_mol, (0, 0))
#deviled = Chem.CanonSmiles(fried)


#ee = Chem.molzip(fried)

    
#combo = Chem.CombineMols(core, r_group1, r_group2, r_group3)
#edcombo.AddBond(5,8,order=Chem.rdchem.BondType.SINGLE)

#FROM HERE - WORK ON ATACHING R GROUPS TO CORES IN EVERY POSSIBLE WAY
#AND COUNTING AMOUNT OF GROUPS
        
#core = Chem.MolFromSmiles(cores[0])
#r_group = Chem.MolFromSmiles(r_groups[0])

#combo = Chem.CombineMols(core,r_group)
#combo
'''
Doesnt join them, just puts them in the same rdkit mol object
'''
#Chem.MolToSmiles(combo)

#edcombo = Chem.EditableMol(combo)
#DrawingOptions.includeAtomNumbers=True
#combo

#edcombo.AddBond(5,8,order=Chem.rdchem.BondType.SINGLE)
'''
adds bond between these 2 specified atoms
'''
#back = edcombo.GetMol()
#back

#core_ring_atoms = get_ring_atoms(core)

def mol_with_atom_index(mol):
    atom_list = []
    for atom in mol.GetAtoms():
        a = atom.SetAtomMapNum(atom.GetIdx())
        b = atom.GetMass()
        print(b)
    return mol

atomsss = mol_with_atom_index(mol)
#### Labels the atoms in the ring 

'''
for b in frag_mol.GetBonds():
    a = (b.GetBeginAtomIdx(),b.GetEndAtomIdx(),
          b.GetBondType(),b.GetStereo())
    #print(a) #- printing this will show type of bond
    b = repr(a)
    print(b)
#This could be useful for fun
''' 

#a = get_ring_atoms(frag_mol)

#m = Chem.MolFromSmiles('C1OC1')
#for atom in m.GetAtoms():
#    print(atom.GetAtomicNum())
    

#c = get_number_of_atoms(frag_mol)


def get_number_atoms(mol):
    smiles = Chem.MolToSmiles(mol)
    print(smiles)
    smiles = smiles.replace("(", "")
    smiles = smiles.replace(")", "")
    smiles = smiles.replace("=", "")
    smiles = smiles.replace("*", "")
    print("egg", smiles)
    no_atoms = len(smiles)
    return(no_atoms)
'''
The above function is probably useless, but it printed out something interesting:
    
    C[CH2:1][CH2:2][OH:3].[cH:4]1[cH:5][cH:6][cH:7][cH:8][cH:9]1
    egg C[CH2:1][CH2:2][OH:3].[cH:4]1[cH:5][cH:6][cH:7][cH:8][cH:9]1


'''

















