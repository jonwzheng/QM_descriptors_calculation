import logging
import os
import json
import time
from rdkit import Chem

REPLACE_LETTER = {"(": "_", ")": "_", "'": "_"}

def mol2xyz(mol):
    return mol.ToXYZ()

def mol2charge(mol):
    return Chem.GetFormalCharge(mol)

def mol2mult(mol):
    num_radical_elec = 0
    for atom in mol.GetAtoms():
        num_radical_elec += atom.GetNumRadicalElectrons()
    return num_radical_elec + 1