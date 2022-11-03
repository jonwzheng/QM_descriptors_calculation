#!/usr/bin/env python
# coding: utf-8

import os
import re
import sys
import shutil

import numpy as np
import pandas as pd
import pickle as pkl

from joblib import Parallel, delayed
from rdkit import Chem

def parser(mol_confs_sdf):

    failed_jobs = dict()
    valid_mol = dict()
    
    mol_id = os.path.basename(mol_confs_sdf).split("_confs.sdf")[0]
    mol_smi = df.loc[df['id'] == mol_id]['smiles'].tolist()[0]
    pre_adj = Chem.GetAdjacencyMatrix(Chem.MolFromSmiles(mol_smi))
    
    failed_jobs[mol_id] = dict()
    valid_mol[mol_id] = dict()
    
    mols = Chem.SDMolSupplier(mol_confs_sdf, removeHs=False, sanitize=True)
    for conf_id, mol in enumerate(mols):
        try:
            xyz = Chem.MolToXYZBlock(mol)
            en = mol.GetProp("ConfEnergies")
            valid_mol[mol_id][conf_id]["ff_xyz"] = xyz
            valid_mol[mol_id][conf_id]["ff_en"] = en
        except:
            failed_jobs[mol_id][conf_id] = "failed"
            
        
    return failed_jobs, valid_mol

input_smiles_path = sys.argv[1]
output_file_name = sys.argv[2]
n_jobs = int(sys.argv[3])

# input_smiles_path = "reactants_products_wb97xd_and_xtb_opted_ts_combo_results_hashed_chart_aug11b.csv"
# n_jobs = 8

df = pd.read_csv(input_smiles_path)
mol_confs_sdf_paths = []
submit_dir = os.getcwd()
for suboutput_folder in os.listdir(os.path.join(submit_dir, "output", "FF_conf", "outputs")):
    for mol_confs_sdf in os.listdir(os.path.join(submit_dir, "output", "FF_conf", "outputs", suboutput_folder)):
        if ".sdf" in mol_confs_sdf:
            mol_confs_sdf_paths.append(os.path.join(submit_dir, "output", "FF_conf", "outputs", suboutput_folder, mol_confs_sdf))

out = Parallel(n_jobs=n_jobs, backend="multiprocessing", verbose=5)(delayed(parser)(mol_confs_sdf) for mol_confs_sdf in mol_confs_sdf_paths)

with open(os.path.join(submit_dir, f'{output_file_name}.pkl'), 'wb') as outfile:
    pkl.dump(out, outfile)

xyz_FF_opt = {}
for failed_dict, success_dict in out:
    if success_dict:
        for mol_id in success_dict:
            xyz_FF_opt[mol_id] = success_dict[mol_id]["ff_xyz"]

with open(f"{output_file_name}_xyz.pkl", "wb") as f:
    pkl.dump(xyz_FF_opt, f)