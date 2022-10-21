# -*- coding: utf-8 -*-
"""
Created on Fri Oct 21 14:12:45 2022

@author: 983045
"""

from functions import *
import smiles_syntax
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem.Draw import MolDrawing, DrawingOptions

DrawingOptions.bondLineWidth=1.8
DrawingOptions.atomLabelFontSize=14
DrawingOptions.includeAtomNumbers=False

mol = Chem.MolFromSmiles("CC(=O)NC1=C(C=C(C=C1)O)O")
#mol = Chem.MolFromSmiles("CC(=O)NC1=C(C=C(C=C1)O)O")
new_mol = cut_bonds_to_core(mol)
print(new_mol)

