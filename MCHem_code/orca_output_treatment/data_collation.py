# -*- coding: utf-8 -*-
"""
Created on Fri Feb 11 18:23:29 2022

@author: Daniel York
"""
import os
import shutil
import glob as glob
import ase as ase
import sys
from ase import Atoms
from ase.build.molecule import molecule
from ase.visualize import view
from ase.io import read,write
from ase.io.xyz import write_xyz
from ase.build import add_adsorbate
from math import pi
from ase.calculators.orca import ORCA
from ase.data.pubchem import pubchem_atoms_search, pubchem_atoms_conformer_search
from ase.collections import g2
from ase.build.molecule import extra
from ase.build import molecule
import numpy as np
import csv
from isomers.functions import *
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import MolFromInchi, MolFromSmiles, MolToSmiles, MolToInchi
from rdkit.Chem import Descriptors


all_files = glob.glob("C:/Users/Daniel York/OneDrive/Documents/mchem_project/isogen1/output_code/refined_seaweed/2-Methyl-2-cyclopenten-1-one/*")
directories = []
molecule_names = []
total_energies = []
homo_energies = []
lumo_energies = []
chemical_hardness_list = []
molecular_weights = []
dipole_moments = []
polarizabilities = []
'''
Empty lists for all the data I want
'''
for everything in all_files:
    isDir = os.path.isdir(everything)
    if isDir == False:
        print(everything, "is not a directory.")
    if isDir == True:
        print("Is this path a directory?", isDir)
        directories.append(everything)
'''
Checks what the files are
Keeps the directories to loop through
'''
for i in range(len(directories)):
    directory = directories[i-1]
    xyz_filepath = directory + "/*.xyz"
    xyz_file = glob.glob(xyz_filepath)
    output_filepath = directory + "/*.out"  
    output_file = glob.glob(output_filepath)  
    
    '''
    Opens required files
    '''
    mol = mol_from_xyz(xyz_file[0])
    smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
    mw = Descriptors.MolWt(mol)
    molecular_weights.append(mw)
    
    
    
    xyz_file = open(xyz_file[0], "r")
    mol = ase.io.read(xyz_file)
    atomic_number_array = mol.numbers
    atomic_num = np.sum(atomic_number_array)
    assert atomic_num % 2 == 0
    homo_num = int(atomic_num/2 - 1)
    lumo_num = homo_num + 1
    print("The homo is:", homo_num, "and the lumo is:", lumo_num)
    xyz_file.close()
    '''
    Finds homo/lumo number
    This is very important for finding the orbital energies!!!!
    '''
    string1 = "TOTAL SCF ENERGY"
    f = open(output_file[0], "r")
    flag = 0
    index = 0
    line_list = []
    for line in f:
        index += 1
        if string1 in line:
            flag = 1
            line_list.append(index)
           # print(index)
            #print(line_list)
    f.close()
    '''
    Finds what lines 'string1' is present in
    Appends these line numbers to a list
    '''
    f = open(output_file[0], "r")
    content = f.readlines()
    lines_to_read = line_list[len(line_list)-1]
    i = lines_to_read
    str_dictionary = repr(content[i+2])
    cut1 = str_dictionary.split("Eh")[1]
    cut2 = cut1.split("eV")[0]
    total_energies.append(cut2)
    '''
    Takes the final line number the string was found
    Cuts up this line to be only the energy
    Appends this to a list of total energies
    '''
    name = content[121]
    new_name = name.split("=")[1]
    molecule_name = new_name.split(".")[0]
    molecule_names.append(molecule_name)
    '''
    Appends molecule name to a list
    '''   
    string2 = "ORBITAL ENERGIES"
    f = open(output_file[0], 'r')
    flag = 0
    index = 0
    line_list = []
    for line in f:
        index += 1
        if string2 in line:
            flag = 1
            line_list.append(index)
            #print(index)
            #print(line_list)
    f.close()
    '''
    Finds what lines 'ORBITAL ENRGIES' is present in
    Appends these line numbers to a list
    '''  
    f = open(output_file[0], 'r')
    content = f.readlines()
    lines_to_read = line_list[len(line_list)-1]
    #print("Lines to read", lines_to_read)
    homo_energy_line = lines_to_read + homo_num + 3
    lumo_energy_line = homo_energy_line + 1
    '''
    Finds the specific lines where the homo/lumo energies are stored 
    Uses the homo/lumo number of the molecule to do this
    '''
    new_line = (content[homo_energy_line])
    new_line_2 = (' '.join(new_line.split()))
    new_line_3 = new_line_2.replace(" ", "#")
    homo_energy = new_line_3.split("#")[3]
    #print("homo energy is:", homo_energy)
    '''
    Finds homo energy value
    The line containing the energy contains random length whitespace
    The random length whitespace is replaced by a single whitespace
    This is replaced with asterisks and used to cut the line to get the energy
    '''    
    new_line = (content[lumo_energy_line])
    new_line_2 = (' '.join(new_line.split()))
    new_line_3 = new_line_2.replace(" ", "#")
    lumo_energy = new_line_3.split("#")[3]
    #print("lumo energy is:", lumo_energy)
    '''
    See above, but instead is for lumo energy
    '''
    lumo_energy = float(lumo_energy)
    homo_energy = float(homo_energy)    
    chemical_hardness_total = (lumo_energy - homo_energy)/2
    rounded_chem_hardness = format(chemical_hardness_total, ".4f")
    #print("Chemical hardness is:", chemical_hardness)
    '''
    Calculates chemical hardness from homo/lumo energies
    '''
    string1 = "Isotropic polarizability"
    file = open(output_file[0], "r")
    for line in file:
        if string1 in line:
            new_line = (' '.join(line.split()))
            polarizability = new_line.split(":")[1]
            polarizabilities.append(polarizability)
    
    string1 = "Magnitude (Debye)"
    file = open(output_file[0], "r")
    for line in file:
        if string1 in line:
            new_line = (' '.join(line.split()))
            dipole_moment = new_line.split(":")[1]
            dipole_moments.append(dipole_moment)
    
    homo_energies.append(homo_energy)
    lumo_energies.append(lumo_energy)
    chemical_hardness_list.append(rounded_chem_hardness)
    '''
    Appends data to lists
    '''
final_dat = []
final_dat.append(molecule_names)
final_dat.append(total_energies)        
final_dat.append(homo_energies)
final_dat.append(lumo_energies)
final_dat.append(chemical_hardness_list)
final_dat.append(molecular_weights)
final_dat.append(dipole_moments)
final_dat.append(polarizabilities)
'''
Creates a 2d array containing all the info for the molecules
'''

with open('2-Methyl-2-cyclopenten-1-one.csv', 'w', encoding='UTF8', newline='') as f:
    writer = csv.writer(f)
    header = ['Molecule name', 'Molecular weight', 'Chemical hardness', 'Dipole moment', 'Polarizability']
    writer.writerow(header)
    for i in range(len(final_dat[0])):
        print("For:", final_dat[0][i-1], "the total energy is:", final_dat[1][i-1], "and homo/lumo energies are:", final_dat[2][i-1], "and", final_dat[3][i-1], "chemical hardness is", final_dat[4][i-1], 'The molecular weight is:', final_dat[5][i-1], "The dipole moment is:", final_dat[6][i-1], "Polarizability", final_dat[7][i-1])
        molecule_name = str(final_dat[0][i-1])
        chemical_hardness = str(final_dat[4][i-1])
        molecular_weight = str(final_dat[5][i-1])
        dipole_moment = str(final_dat[6][i-1])
        polarizability = str(final_dat[7][i-1])
        data = ['%s' % molecule_name, '%s' % molecular_weight, '%s' % chemical_hardness, '%s' % dipole_moment, '%s' % polarizability]
        writer.writerow(data)



   
   
  
      
     
     