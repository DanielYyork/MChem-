# -*- coding: utf-8 -*-
"""
Created on Tue Mar 22 18:41:30 2022

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

def isRingAromatic(mol, bondRing):
        for id in bondRing:
            if not mol.GetBondWithIdx(id).GetIsAromatic():
                return False
        return True
    
ring_sizes = []
for i in range(3, 8):
    ring_size = '[r' + str(i) + ']'
    ring_sizes.append(ring_size)
    
atoms = ["C", "c", "O", "o", "N", "n"]
functionals = ["[O;H1]", "[O;D1]", "[O;D2]", "[CX3]=[OX1]"]

'''
Creates a list of things to search for - in this case it is differently sized rings
C = aliphatic carbon
c = aromatic carbon
functionals = [simple hydroxy oxygen, 1-connected (hydroxy or hydroxide) oxygen, 2-connected (etheric) oxygen, carbonyl group]
'''
name = "2-hydroxy-3-methyl-2-cyclopenten-1-one"
original_smiles = "CC1=C(C(=O)CC1)O"
structures_to_search_for = []

mol = Chem.MolFromSmiles(original_smiles)
for ring_size in ring_sizes:
    pattern_smarts = Chem.MolFromSmarts(ring_size)
    if mol.HasSubstructMatch(pattern_smarts) == True:
        print("The initial smilesstring is:", original_smiles, " and it has a", ring_size, "ring")
        structures_to_search_for.append(ring_size)
ri = mol.GetRingInfo()
number_of_rings = len(ri.AtomRings())

atoms_in_original_smiles = []
for atom in atoms:
    pattern_smarts = Chem.MolFromSmarts(atom)
    if mol.HasSubstructMatch(pattern_smarts) == True:
        print("aluaha snahckbar")
        atoms_in_original_smiles.append(atom)
structures_to_search_for.append(atoms_in_original_smiles)

oxygen_functionals_in_original_smiles = []
for functional in functionals:
    pattern_smarts = Chem.MolFromSmarts(functional)
    if mol.HasSubstructMatch(pattern_smarts) == True:
        print("Im already tracer")
        oxygen_functionals_in_original_smiles.append(functional)
structures_to_search_for.append(oxygen_functionals_in_original_smiles)

sp = (sum((x.GetHybridization() == Chem.HybridizationType.SP) for x in mol.GetAtoms()))
sp2 = (sum((x.GetHybridization() == Chem.HybridizationType.SP2) for x in mol.GetAtoms()))
sp3 = (sum((x.GetHybridization() == Chem.HybridizationType.SP3) for x in mol.GetAtoms()))
structures_to_search_for.append(sp)
structures_to_search_for.append(sp2)
structures_to_search_for.append(sp3)
'''
Finds patterns in the original molecule
Format of "structures_to_search_for":
    [ring size, sp#, sp2#, sp3#, aromaticity]
Finds the ring size and hybridaization (i.e. sp indicates a triple bond)
'''
aromatic_rings = 0
for i in range(len(list(ri.AtomRings()))):
    if isRingAromatic(mol, ri.BondRings()[i-1]) == True:
        aromatic_rings += 1        
        print(original_smiles, "has aromatic ring")
        structures_to_search_for.append("Aromatic ring")
    else:
        print(original_smiles, "doesnt have aromatic ring")
        structures_to_search_for.append("Non aromatic ring")
'''
Iterates over each ring and checks to see if these rings are aromatic
'''
'''
Below is where is begins
'''       
molecules_with_rings = [[],[],[],[],[]]

file = open('C6H8O2.smi', 'r')

line_list = []
for line in file:
    stripped_line = line.strip()
    line_list.append(stripped_line)
'''
Appends each smilesstring made by the isoemr generator to a list
molecules_with_rings = [[3],[4],[5],[6]]
'''   
for i in range(len(line_list)-1):
    mol = Chem.MolFromSmiles(line_list[i])
    for ring_size in ring_sizes:
        pattern_smarts = Chem.MolFromSmarts(ring_size)
        if mol.HasSubstructMatch(pattern_smarts) == True:
            #print(line_list[i], "has a", ring_size, "ring")
            place_to_put_in_big_list = ring_size[2]
            molecules_with_rings[int(place_to_put_in_big_list)-3].append(line_list[i])
'''
Line_list is every smiles string relating to each isomer created
For each of these, a check for what rings it contains is done
The -3 is there simply because index 0 corresponsds to 3 membered rings as 1 and 2 memeber rings are impossible
'''
print("There are:", len(line_list), "isomers for C6H8O2")
print("There are:", len(molecules_with_rings[0]), "with a 3-membered ring")
print("There are:", len(molecules_with_rings[1]), "with a 4-membered ring")
print("There are:", len(molecules_with_rings[2]), "with a 5-membered ring")
print("There are:", len(molecules_with_rings[3]), "with a 6-membered ring")
# NEED TO GET IT TO CALCULATE HOW MANY RINGS PER ISOMER?
# HOW MANY HAVE DIFFERENT AMOUNTS OF RINGS (e.g, cna have 3 and 6 membered rings)
'''
The section below takes the identifed ring and explores that list of that ring size
'''
refined_isomer_list = [] #

if number_of_rings == 1:
    print("Initial smilesstring contains 1 ring, which is", structures_to_search_for)
    name = repr(structures_to_search_for[0])
    new_name = name.split('r')[1]
    final_ring_size = int(new_name.split(']')[0])
    molecules_to_refine = molecules_with_rings[final_ring_size-3] #This line selects the isomers with the same ring size as the initial molecule
    '''
    Searches through the list of 5 memebered ring isomers
    '''
a = 0      
for i in range(len(molecules_to_refine)):    
    mol = Chem.MolFromSmiles(molecules_to_refine[i-1])
    ri = mol.GetRingInfo()
    number_of_rings_in_isomers = len(ri.AtomRings())
    number_of_molecules_removed = 0
    if number_of_rings_in_isomers > number_of_rings:
        a += 1
    else:
        refined_isomer_list.append(molecules_to_refine[i-1])
'''
Sorts though five membered rings and pulls out any isomers with a different number of rings to the original molecule
'''
refined_isomer_list_2 = []
for i in range(len(refined_isomer_list)):
    mol = Chem.MolFromSmiles(refined_isomer_list[i-1])
    sp = (sum((x.GetHybridization() == Chem.HybridizationType.SP) for x in mol.GetAtoms()))
    sp2 = (sum((x.GetHybridization() == Chem.HybridizationType.SP2) for x in mol.GetAtoms()))
    sp3 = (sum((x.GetHybridization() == Chem.HybridizationType.SP3) for x in mol.GetAtoms()))
    if structures_to_search_for[4] != 0 and sp == 0:
        continue
    if structures_to_search_for[5] != 0 and sp2 == 0:
        continue
    refined_isomer_list_2.append(refined_isomer_list[i-1])   
        
refined_isomer_list_2 = set(refined_isomer_list_2)
'''
Compares hybridization of isomers to that of the initial molecule
Quick way to discern sigma, pi, triple bond nature
If initial compound has triple bond (or double)
Isomers with that bond type are kept
'''
refined_isomer_list_3 = []
refined_isomer_list_2 = list(refined_isomer_list_2)
for i in range(len(refined_isomer_list_2)):
    mol = Chem.MolFromSmiles(refined_isomer_list_2[i-1])
    original_atoms = structures_to_search_for[1]
    for atom in original_atoms:
        pattern_smarts = Chem.MolFromSmarts(atom)
        if mol.HasSubstructMatch(pattern_smarts) == True:
            refined_isomer_list_3.append(refined_isomer_list_2[i-1])
refined_isomer_list_3 = set(refined_isomer_list_3)
          
refined_isomer_list_4 = []
refined_isomer_list_3 = list(refined_isomer_list_3)
for i in range(len(refined_isomer_list_3)):
    mol = Chem.MolFromSmiles((refined_isomer_list_3[i-1]))
    original_functionals = structures_to_search_for[2]
    for functional in original_functionals:
        pattern_smarts = Chem.MolFromSmarts(functional)
        if mol.HasSubstructMatch(pattern_smarts) == True:
            refined_isomer_list_4.append(refined_isomer_list[i-1])
refined_isomer_list_4 = set(refined_isomer_list_4)
refined_isomer_list_4 = list(refined_isomer_list_4)
    
structures_to_search_for_refined = structures_to_search_for[2]

def functional_match(initial_smiles, list_orig_ox, list_iso):   
   
    
    new_list = []
    for i in range(len(list_iso)):
       
        mol = Chem.MolFromSmiles(list_iso[i-1])
        a = len(list_orig_ox)
        b = 0
        for thing in list_orig_ox:
            pattern_smarts = Chem.MolFromSmarts(thing)
            if mol.HasSubstructMatch(pattern_smarts) == True:
                a += 1
        if a == b:
            new_list.append(list_iso[i-1])
    return(new_list)
            
print(functional_match("CC1=C(C(=O)CC1)O", oxygen_functionals_in_original_smiles, refined_isomer_list_4))
      
        
    
        


'''
Initial smiles = smilestring of original molecule
list = list of strucures to search for


        
    
    

            
            
         
print(len(refined_isomer_list))


    
    for i in range(len(list(ri.AtomRings()))):
        if isRingAromatic(mol, ri.BondRings()[i-1]) == True:
            aromatic_rings += 1        
            print(molecules_to_refine[i-1], "has aromatic ring")
        #else:
            #print(molecules_to_refine[i-1], "hasmt aromatic ring")
     
    #else:
     #   refined_isomer_list.append(molecules_to_refine[i-1])
print("There number of molecules with same ring size and number of rings as initial smiles is:", len(refined_isomer_list))









Goes through an amount of specified molecules and outputs what size ring each one has
Appends to 2d array
with [0][*] = 3 membered rings
with [1][*] = 4 membered rings
with [2][*] = 5 membered rings
with [3][*] = 6 membered rings
with [4][*] = 7 membered rings

NEED TO:
- MAKE IT IDENTIFY WHAT SIZE RINGS ARE IN ORIGINAL MOLECULE
- OR LENGTH OF ALIPHATIC CHAIN
- PULL OUT ISOMERS THAT SHARE THIS CHARECTERISTIC

Line_list has every smilesstring created by the isomer generator
 
    
print(line_list[:11])

molecules_with_5rings = []
molecules_with_6rings = []
molecules_with_aromatic_rings = []
a = 0  
pattern_smarts = Chem.MolFromSmarts('[r5]')
for molecule in line_list:
    mol = Chem.MolFromSmiles(molecule)
    if mol.HasSubstructMatch(pattern_smarts):
        molecules_with_rings.append(molecule)
        ri = mol.GetRingInfo()
        aromatic_rings = 0
        for i in range(len(list(ri.AtomRings()))):
            if isRingAromatic(mol, ri.BondRings()[i-1]) == True:
                aromatic_rings += 1
                a += 1
                molecules_with_aromatic_rings.append(molecule)
               

            

# You can interrogate the RingInfo object to tell you the atoms that make up each ring:
print(ri.AtomRings())

# or the bonds that make up each ring:
print(ri.BondRings())

# To detect aromatic rings, I would loop over the bonds in each ring and
# flag the ring as aromatic if all bonds are aromatic:







pattern = Chem.MolFromSmiles('N')
# Checks if it has nitrogen
pattern_smarts = Chem.MolFromSmarts('[r]')
# Checks if it has a ring 'r5' means 5-membered ring and so on...

for mol in mol_list:
    print(mol.HasSubstructMatch(pattern))
    
for mol in mol_list:
    print(mol.HasSubstructMatch(pattern_smarts))
    
    '''