# -*- coding: utf-8 -*-
"""
Created on Tue Oct 18 11:11:54 2022

@author: 983045
"""
from rdkit import Chem
test_mol = Chem.MolFromSmiles("CC(=O)NC1=C(C=C(C=C1)O)O")

def get_bonds(mol):
    all_bonds = []
    for b in mol.GetBonds():
        a = (b.GetBeginAtomIdx(),b.GetEndAtomIdx(),
              b.GetBondType(),b.GetStereo())
        #print(a) #- printing this will show type of bond
        b = repr(a)
        bond_atoms = []
        atom_one = a[0]
        atom_two = a[1]
        bond_type = repr(a[2]).split('.')[4]
        bond_atoms.append(atom_one)
        bond_atoms.append(atom_two)
        bond_atoms.append(bond_type)
        all_bonds.append(bond_atoms)
    return(all_bonds)

a = get_bonds(test_mol)
'''
The 'get_bonds' function takes an RDkit mol object and returns a 2d array
detailing each bond. It provides the numbers of atoms either end of the bond 
and the type of bond.
'''
def get_ring_atoms(mol):
    ri = mol.GetRingInfo()
    c = ri.AtomRings()
    no_of_rings = len(c)
    ring_atoms = set()
    for i in range(no_of_rings):
        ring_atom = c[i]
        for j in range(len(ring_atom)):
            ring_atoms.add(c[i][j])
    return(ring_atoms)

#b = get_ring_atoms(test_mol)
'''
The 'get_ring_atoms' function takes an RDkit mol object and returns a set. This
set contains the numbers of the atoms that comprise rings within the molecule.

Returns a set of which atoms make up rings within the molecule.
'''
def fragment_simple(mol, atom1, atom2):
    rwmol = Chem.RWMol(mol)
    rwmol.RemoveBond(atom1, atom2)
    wildcard1 = rwmol.AddAtom(Chem.Atom(0))
    wildcard2 = rwmol.AddAtom(Chem.Atom(0))
    rwmol.AddBond(atom1, wildcard1, Chem.BondType.SINGLE) 
    rwmol.AddBond(atom2, wildcard2, Chem.BondType.SINGLE) 
    return rwmol.GetMol()
'''
The 'fragment_simple' takes an RDkit mol object and the atom numbers of atoms
either side of the bond. The bond between these two atoms is then split and wildcard
atoms are added to each end of the cut bond - which are useful for reforming the molecule.

Returns a fragmented mol

This function is taken from:
    http://www.dalkescientific.com/writings/diary/archive/2016/08/14/fragment_chiral_molecules.html

NOTE: This only works properly for achiral molecules and chirality should be preserved
somehow - allowing for the construction of even more unique molecules.
'''
def cut_bonds_to_core(mol):
    a = get_bonds(mol)
    b = get_ring_atoms(mol)
    bond_pairs = []   
    for i in range(len(a)):
        atom_pair = []
        atom_one = a[i][0]
        atom_two = a[i][1]
        if atom_one in b:
            if atom_two in b:
                #print("atoms makeing ring bond", atom_one, atom_two)
                pass
            elif atom_two not in b:
                #print("ring atom to non ring atom bond", atom_one, atom_two)
                atom_pair.append(atom_one)
                atom_pair.append(atom_two)
                bond_pairs.append(atom_pair)
        if atom_two in b:
            if atom_one in b:
                #print("skipping")  
                pass
            elif atom_one not in b:
                #print("ring atom to non ring atom bond", atom_one, atom_two )
                atom_pair.append(atom_one)
                atom_pair.append(atom_two)
                bond_pairs.append(atom_pair) # This first part of the function identifies which bonds to be cut                
    for i in range(len(bond_pairs)):
        wildcard_replacer = 80
        wildcard = '*'
        wildcard_branch = '(*)'
        if i == 0:
            frag_time = fragment_simple(mol, bond_pairs[i][0], bond_pairs[i][1])
            print("this is frag time", Chem.MolToSmiles(frag_time), type(Chem.MolToSmiles(frag_time)))
            #if wildcard in frag_time:
             #   frag_time.replace(wildcard, str(wildcard_replacer))
        else:
            frag_time = fragment_simple(frag_time, bond_pairs[i][0], bond_pairs[i][1])
            print("this is frag time", Chem.MolToSmiles(frag_time))
    cut_mol = Chem.MolToSmiles(frag_time) 
    return(cut_mol) 
'''
The 'cut_bonds_to_core' function takes a single RDkit mol object and returns a SMILES
where any bonds to a ring core have been cut. This will only work with ring core structures
and aliphatic molecules with chains must be considered differently. 

Returns a smiles string of the snipped up molecule.

NOTE: No function for removing branches from aliphatic molecules currently exists.
'''
def has_rings(mol):
    ring_functionals = ["[r5]", "[r6]", "[r7]", "[r8]"]
    ring_names = ["5-membered ring", "6-membered ring", "7-membered ring", "8-membered ring"]
    ring_groups = []
    ring_groups.append(ring_names)
    ring_groups.append(ring_functionals)
    ring_groups_in_mol = []
    for i in range(len(ring_groups[0])):
        pattern_smarts = Chem.MolFromSmarts(ring_groups[1][i])
        if mol.HasSubstructMatch(pattern_smarts) == True:
            ring_groups_in_mol.append(ring_groups[0][i])
    if len(ring_groups) != 0:
        return(ring_groups_in_mol)
    if len(ring_groups) == 0:
        no_list = ["N"]
        return(no_list)


def get_atom_set(mol):
    counter = 0
    for atom in mol.GetAtoms():
        counter += 1
    atom_set= set()
    for i in range(counter):
        atom_set.add(i+1)
    return(atom_set)

    
