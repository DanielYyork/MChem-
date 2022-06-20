# -*- coding: utf-8 -*-
"""
Created on Thu Apr 21 10:06:31 2022

@author: Daniel York
"""
from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Draw
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem
from rdkit import DataStructs
import numpy as np
from PIL import Image
from rdkit.Chem import rdmolops
import os
import glob as glob

def unique(listt):
    new_set = set(listt)
    new_list = list(new_set)
    return(new_list)
'''
Unique
ARGUMENTS: (list)
Returns the list, but only with unique elements
i.e. removes any duplicates
'''
name = "4-aminophenol"
initial_smiles = "C1=CC(=CC=C1N)O"


heterocycle_functionals = ["[oX2r5]","[nX2r5]", "[sX2r5]", "[oX2r6]","[nX2r6]", "[sX2r6]", "[n]"]
carbonyl_functionals = ["[NX3][CX3](=[OX1])[#6]", "[CX3H1](=O)[#6]", "[#6][CX3](=O)[#6]", "[CX3](=O)[OX2H1]", "[OD2]([#6])[#6]", "[#6][CX3](=O)[OX2HO][#6]", "[#6][CX3](=O)[#6].[OD2]([#6])[#6]", "[#0][CX3](=O)[#0]"]
#amine_functionals = ["NX3;H2,H1;!$(NC=[!#6])", "NX3;H1;!$(NC=[!#6])", "[NX3][$(C=C,$(cc)"]
nitrogenous_functionals = ["[Nx1]#[CX2]", "[$([NX3](=O)=O),$([NX3+](=O)[O-])][!#8]", "[NX2]=[OX1]", "[NX3&H2]", "[NX3&H1]"]
hydroxyl_groups = ["[#6][OX2H]", "[OX2H][CX3]=[OX1]", "[OX2H][cX3]:[c]"]
#alkane_chains = ["[RO;D2][RO;D2]", "[RO;D2][RO;D2][RO;D2]", "[RO;D2][RO;D2][RO;D2]"]
alkane_chains = ["[CH2]", "[CH2][CH2]", "[CH2][CH2][CH2]", "[CH2][CH2][CH2][CH2]", "[CH2][CH2][CH2][CH2][CH2]", "[CH2][CH2][CH2][CH2][CH2][CH2]", "[CH2][CH2][CH2][CH2][CH2][CH2][CH2]", "[CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2]", "[CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2]"]
carbonness = ["[CD2]", "[CD3]", "[CD4]", "[CH3]", "[CX3&H0]"]
#any_chains = ["[RO;D2]~[RO;D2]", "[RO;D2]~[RO;D2]~[RO;D2]", "[RO;D2]~[RO;D2]~[RO;D2]"]
ring_sizes = ['[r3]', '[r4]', '[r5]', '[r6]', '[r7]'] 
functional_lib = []
functional_lib.append(heterocycle_functionals)
functional_lib.append(carbonyl_functionals)
#functional_lib.append(amine_functionals)
functional_lib.append(nitrogenous_functionals)
functional_lib.append(hydroxyl_groups)
functional_lib.append(carbonness)
functional_lib.append(alkane_chains)
#functional_lib.append(any_chains)
functional_lib.append(ring_sizes)
'''
Above code is creation of a library of functionals to search for
- Any functionals can be specified here, as long as they follow SMARTS rules
heterocycle functionals = [furan, pyrrole, etc..., n (in ring?!)]
carbonyls = [Amide, Aldehyde, Ketone, cooh, ether, ester, carbonyl, dicarbonyl (will search for two carbonyls connected by a carbon atom, but either any atoms can be the r groups)]
amines = [Primary amine, secondary amine, enamine, primary amine, secondary amine]
nitrogenous = [Nitrile, nitro, nitroso]
hydroxyl = [in alcohol, in cooh, phenol]
carboness = [carbon bonded to 2 carbon, .. to 3 carbon, .. to 4 carbon]
alkane chains = [ethyl, propyl, butyl] - ONLY searches alkanes
any chains = [ethyl, propyl, butyl] -- Searches, alkanes, alkenes and alkynes
ring_ sizes = [3, 4, 5, 6, 7]
NOTE: hashed out groups may not work currently
'''
def functional_search(functionals, mol):
    # Functionals = list of functionals/things to search for
    # Molecule = rdkit mol object
    init_funct = []
    for functional in functionals:
        pattern_smarts = Chem.MolFromSmarts(functional)
        if mol.HasSubstructMatch(pattern_smarts) == True:
            init_funct.append(functional)
    return(init_funct)
'''
functional_search:
Arguments: (list of functionals, rdkit mol object)
This is used inside the larger struct search
The list of functionals are pulled from a created functional library
Returns the functionals that are in mol object (if the functional is in the functional library)
'''
def get_functional_amount(listt, mol):
    init_funct_amount = []
    for thing in listt:
        search = Chem.MolFromSmarts(thing)
        funct_amount = (len(mol.GetSubstructMatches(search)))        
        init_funct_amount.append(funct_amount)
    return(init_funct_amount)
'''
get_functional_amount
Arguments: (list of funtionals #from 'functional_search#', rdkit mol object)
Returns the amount of each identified functional in the original molecule
'''
def struct_search(smiles, functional_lib):
    mol = Chem.MolFromSmiles(smiles)
    '''
    This will return a list with key molecular structures and patterns
    structures_to_search_for = [[ring size], [atoms], [double bond], [triple bond] [oxygen functionals], [nitrogen functionals]]
    '''
    struct_list = []
    struct_amount = []
    for i in range(len(functional_lib)):
        functional_list = functional_lib[i-1]
        init_funct = functional_search(functional_list, mol)
        struct_list.append(init_funct)
        init_funct_amount = get_functional_amount(init_funct, mol)
        struct_amount.append(init_funct_amount)
    '''
    This loop goes through every functional in the functional library;
    and adds it to struct_list and struct_amount 
    (i.e. finds functionals and how many of each functional - obvilosuly it can only find functional in the library)
    '''
    final_list = []
    final_list.append(struct_list)
    final_list.append(struct_amount)
    '''
    Final_list contains what groups/functionals in the original molecule as well as the amount of that functional
    Structure:
        Final_list = [[functionals], [functional amounts]]
        Functionals = [[ring size], [atoms], [carbonyls]]
        Functional amounts = [[amount of each ring], [amount of each atom], [amount of each carbonyl]]
    '''
    return(final_list)
'''
Struct_search
Arguments: (smilestring, list of functionals to search for)
The smilestring is the molecule you want have things identified in
The list of functionals to search for must be created, but this can contain any SMARTS pattern
'''
def isomer_refine_by_functional(isomer_list, single_functional, functional_amount):
# Refines a list of isomers by a comparing to a list of required functionals
# The isomer list is the lsit of isomers genereated by MAYGEN
# The funct _ident_list is a part of the funct_ident array   (i.e. a list of a type of fuinctional icdentified in the initial molecule)     
    refined_list = []   
    for i in range(len(isomer_list)):
        mol = Chem.MolFromSmiles(isomer_list[i-1])
        pattern_smarts = Chem.MolFromSmarts(single_functional)
        if mol.HasSubstructMatch(pattern_smarts) == True:
            search = Chem.MolFromSmarts(single_functional)
            funct_amount = (len(mol.GetSubstructMatches(search))) 
            if funct_amount == functional_amount:
                refined_list.append(isomer_list[i-1])
    return(refined_list)
funct_ident = struct_search(initial_smiles, functional_lib)
'''
Isomer_refine_by_functional
ARGUMENTS: List of isomers/molecule, A single functional (SMARTS string), integer (amount of functional want identified)
This takes a functional with a specified amount of that functional;
and filters out any molecules/isomers from a list with those charecteristics.
This can be adapted to remove the amount and just search by functional OR macrostrutures
'''
isomers = []
#file = open("C:/Users/Daniel York/OneDrive/Documents/mchem_project/isogen1/isomer_treatment/2-hydroxy-3-methyl-2-cyclopenten-1-one/C6H8O.smi", "r")
file = open("C:/Users/Daniel York/OneDrive/Documents/mchem_project/isogen1/isomer_treatment/2-hydroxy-3-methyl-2-cyclopenten-1-one/C6H7NO.smi")
for line in file:
    stripped_line = line.strip()
    isomers.append(stripped_line)

all_store = [ [] for _ in range(len(funct_ident[0])) ] 
print("Initial number of isomers is:", len(isomers))
for i in range(len(funct_ident[0])):    
    if i == 0:
        refine_store = isomers
    for j in range(len(funct_ident[0][i-1])): # functional group
        functional = funct_ident[0][i-1][j-1]
        print(functional)
        functional_amount = funct_ident[1][i-1][j-1] 
        refine = isomer_refine_by_functional(refine_store, functional, functional_amount)
        refine_store = refine
        print(len(refine_store), "is the number of isomers after filter:", i+1, ".", j+1)
        if i != 0:
            print(i)
            all_store[i].append(refine_store)
             
refine = unique(refine)
print("Final number of unique isomers is:", len(refine))

with open('refined_isomers.smi', 'w') as f:
    for item in refine:
        f.write("%s\n" % item)
    f.close()
        
with open('refined_isomers.names', 'w') as f:
    for i in range(len(refine)):
        new_name = name + str(i+1)
        f.write("%s\n" % new_name)
    f.close()
        
        
    
