# -*- coding: utf-8 -*-
"""
Created on Thu Feb 10 12:59:38 2022

@author: Daniel York

"""
import shutil
import glob as glob
import os
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

def append_new_line(file_name, text_to_append):
    """Append given text as a new line at the end of file"""
    # Open the file in append & read mode ('a+')
    with open(file_name, "a+") as file_object:
        # Move read cursor to the start of file.
        file_object.seek(0)
        # If file is not empty then append '\n'
        data = file_object.read(100)
        if len(data) > 0:
            file_object.write("\n")
        # Append text at the end of file
        file_object.write(text_to_append)
# ^^ This part of code need to be taken out and called - instead f being here

#smilefiles = list(glob("C:/Users/Daniel York/OneDrive/Documents/mchem_project/isogen1/2_feb/smiles/*.smi"))
#namefiles = list(glob("C:/Users/Daniel York/OneDrive/Documents/mchem_project/isogen1/2_feb/molecule_names/*.names"))
# I want to use 'glob' as above but it's not working properly 

smilefiles = glob.glob("*.smi")
namefiles = glob.glob("*.names")
# Gathers relevant files - currently a hotfix - see above
total_files = len(smilefiles)
file_number = 0
# This is variable that gets updated every time the loop runs - tells the code which file to open
for i in range(total_files):
    if i == total_files:
        break
    namefile = namefiles[file_number]
    file = open(namefile)
    name_list = []
    # Creates an empty list
    for line in file:
        stripped_line = line.strip()
        line_list = stripped_line.split()
        name_list.append(line_list)
    # Puts all molcule names in a list
    smilefile = smilefiles[file_number]
    file = open(smilefile)
    name = smilefile
    str_dictionary = repr(name) 
    file_to_make = name.split(".")[0] # Takes alage.smi and cuts the string down to 'algae' - these will be names of the directories             
    # Update this path to folder you want to save in
    parent_path = 'C:/Users/Daniel York/OneDrive/Documents/mchem_project/isogen1/input_creation' # This is the parent directory where the folders containing input files will be created
    directory = file_to_make                      
    file_number += 1
    namecount = 0 # Will tell the next loop which molecule name to pull from 'name_list' - required to name the input file
    for line in file: # This goes through each smiles string in the folder containing each smiles string (e.g. algae.smi)
        stripped_line = line.strip()
        str_dictionary = repr(stripped_line)  
        print(stripped_line)
        new_file = open("%s.missingfile" % file_to_make, "a")
        try:            
            name = pubchem_atoms_search(smiles=stripped_line)  
            # Searches for the smilestring in pubchem
            # Will create an xyz file for that smilestring if it could be found
        except ValueError:
            print("Could not find compound for: %s" % str_dictionary)
            append_new_line("%s.missingfile" % file_to_make, "%s" %str_dictionary) 
            str_dictionary = repr(name_list[namecount])
            append_new_line("%s.missingfile" % file_to_make, "%s" %str_dictionary)
            # Appends any molecules not able to be found by pubchem to a seperate file - inputs must be made manually for these
            # Creates this file even if all smiles are found - not neccessary
            #saves this file to the same directory of all of the inputfiles
            continue
        molecule = repr(name_list[namecount])
        molecule2 = molecule.split("'")[1]
        # Snips up the name in the list into just molecule name - there is definitely a better way to do it, just a hotfix
        print(name)
        write('%s.xyz' % molecule2, name)
        namecount += 1
        os.path.join(parent_path, directory)
        if not os.path.isdir(directory): # Makes the algae/seaweed/ect.. dir
            os.mkdir(directory)
        xyz_parentdir = parent_path + "/" + directory
        xyz_directory = xyz_parentdir + "/XYZ"
        os.path.join(xyz_parentdir, xyz_directory)
        print
        if not os.path.isdir(xyz_directory):
            os.mkdir(xyz_directory)
        filename = "%s.xyz" % molecule2
        newpath = shutil.copy(filename, xyz_directory)
                
    filename = "%s.missingfile" % file_to_make
    missingfile_path = os.path.join(parent_path, filename)    
    if not os.path.isdir(directory):
        os.mkdir(directory)            
    newpath = shutil.copy(filename, directory)
        # Copies the file containing any missing smile string to the relevant directory (e.g. algae.missingfile will be in the algae directory - althoguh these files need to be removed form the parent dierectory after this)            
        
    xyzfiles = glob.glob('*.xyz')
    '''
    for xyzfile in xyzfiles:
        f = open(xyzfile, "r")
        f.close()
        os.remove(xyzfile)
    
    

    for xyzfile in xyzfiles:
       file = open(xyzfile, "r")
        line_count = 0
        for line in file:
            if line != "\n":
                line_count += 1
        file.close()
            
        
    # Counts lines of all xyz files
'''
    for xyzfile in xyzfiles:
        mol = ase.io.read(xyzfile)
        atomic_number_array = mol.numbers
        atomic_num = np.sum(atomic_number_array)
        assert atomic_num % 2 == 0
        homo_num = int(atomic_num/2 - 1)
        lumo_num = homo_num + 1
        homo_num 
        name = xyzfile.split(".")[0]
        print("For %s" % name)
        print("The homo is", homo_num)
        print("The lumo is", lumo_num) 
        orbitals = str(homo_num) + "." + str(lumo_num)
        f = open("%s.orbitals" % name, "w")
        f.writelines(orbitals)
        f.close()
        # Finds the homo and lumo numbers
        
    for xyzfile in xyzfiles:
        text = open(xyzfile).readlines()
        f = open(xyzfile, "w")
        f.writelines(text[2:])
        f.close()       

    for xyzfile in xyzfiles:
        f = open(xyzfile, "r")
        name = repr(xyzfile)
        middle_name = name.split(".")[0]
        final_name = middle_name.split("'")[1]
        g = open("%s.fixedxyz" % final_name, "w")
        for line in f:
            if line.strip():
                g.write("\t".join(line.split()[:4]) + "\n")
        g.close()            
        f.close() 
     
    fixed_files = glob.glob("*.fixedxyz")        
    # Removes first two lines from xyz file - the first 2 lines are not coordinates       

    for xyzfile in fixed_files:
        #lines = ['! B3LYP def2-TZVP D3BJ keepdens opt', '%scf', 'MaxIter 1000', 'end', '%output', 'Print[P_Hirshfield] 1', 'end', '%elprop', 'Polar 1', 'end', '%plots', 'dim1 100', 'dim2 100', 'dim3 100', 'Format Gaussian_Cube', 'ElDens("filename.dens.cube");', 'end', '%pal', 'nprocs 10', 'end', '*xyz 0 1']
        name = xyzfile
        str_dictionary = repr(name)
        new_name = name.split(".")[0] # Takes xyzfile name as a string and cuts of '.xyz'       
        f = open("%s.inp" % new_name, "a") # Creates an input file with name being isomer smiles 
        f.write
        f.close()
       
    # Creates a text file (input file) for each xyzfile - this is 'naked' and without coordinated at this point  
    

    inputfiles = glob.glob('*.inp')    

    for inputfile in inputfiles:
        
        name = inputfile
        str_dictionary = repr(name)
        new_name = name.split(".")[0]  
        f = open("%s" % name, "a")
        g = open("%s.orbitals" % new_name, "r")
        orbitals = g.readlines()[0]
        homo_num = orbitals.split(".")[0]
        lumo_num = orbitals.split(".")[1]
        lines = ['! B3LYP def2-tzvp D3BJ keepdens opt', '%scf', ' MaxIter 1000', 'end', '%output', ' Print[ P_Hirshfeld] 1', 'end', '%elprop', ' Polar 1', 'end', '%plots', ' dim1 100', ' dim2 100', ' dim3 100', ' Format Gaussian_Cube', ' ElDens("%s.dens.cube");' % new_name, ' MO("homo.cube", %s, 0);' % homo_num, ' MO("lumo.cube", %s, 0);' % lumo_num, 'end', '%pal', ' nprocs 10', 'end', '*xyz 0 1']            
        f.write('\n'.join(lines))  # writes a file with each specified string above being written to a new line   
        #f.close() 
        g.close()
     
    # Writes the 'skeleton' input file with each file made in the last loop    
    # Change in line 138 to change whats in the input
    
    for xyzfile in fixed_files:
        file = open(xyzfile, "r")
        name = xyzfile
        str_dictionary = repr(name)
        new_name = name.split(".")[0]  
        for line in file:
                stripped_line = line.strip()       
                append_new_line("%s.inp" % new_name, stripped_line)          
                #f.close()
        file.close()
                
    # Opens each xyz file and appends coordinates (from corresponding xyz file) to correct file            
                
    for inputfile in inputfiles:
        f = open(inputfile, "a")
        end = "*"
        append_new_line(inputfile, end)        
        f.close()

    # Opens each input file and appends final asterisk 
       
    for inputfile in inputfiles:
        filename = inputfile
        file_path = os.path.join(parent_path, filename)
        if not os.path.isdir(directory):
            os.mkdir(directory)            
        newpath = shutil.copy(filename, directory)
        #newpath = shutil.copy(file_to_make, directory)
        
        
    # Makes a directory for each molecule and puts it isomer input files in there    
        
   
    # The following lines of code clean the parent directory of input and xyz files (input files now in folders and xyz files no longer required)
    # Currently not working for 'missingfiles'
    inputfiles = glob.glob("*.inp") 
    orbital_files = glob.glob("*.orbitals")
    '''
    for orbital_file in orbital_files:
        file = open(orbital_file)
        file.close()
        os.remove(orbital_file)
    '''
    for inputfile in inputfiles:
        file = open(inputfile)
        file.close()
        os.remove(inputfile) 
        
    xyzfiles = glob.glob("*.xyz") 

    for xyzfile in xyzfiles:
        file = open(xyzfile)
        file.close()
        os.remove(xyzfile)
        
    for xyzfile in fixed_files:
        file = open(xyzfile)
        file.close()
        os.remove(xyzfile)
    
    orbital_files = glob.glob("*.orbitals")
    for orbital_file in orbital_files:
        os.remove(orbital_file)
    #missingfiles = glob.glob("*missingfile") 

    #for missingfile in missingfiles:
     #   file = open(missingfile)
      #  file.close()
       # #os.remove(missingfile)
 

    
    
    
    
    
    