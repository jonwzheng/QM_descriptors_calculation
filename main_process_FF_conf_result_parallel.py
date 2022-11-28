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
from rdmc.mol import RDKitMol

def parser(mol_id):

    ids = str(int(int(mol_id.split("id")[1])/1000)) 
    mol_confs_sdf = os.path.join(submit_dir, "output", "FF_conf", "outputs", f"outputs_{ids}", f"{mol_id}_confs.sdf")
    if os.path.isfile(mol_confs_sdf):
        failed_jobs = dict()
        valid_mol = dict()

        mol_smi = mol_id_to_smi[mol_id]
        pre_mol = RDKitMol.FromSmiles(mol_smi)
        pre_adj = pre_mol.GetAdjacencyMatrix()
        
        failed_jobs[mol_id] = dict()
        valid_mol[mol_id] = dict()
        
        mols = RDKitMol.FromFile(mol_confs_sdf)
        for conf_id, mol in enumerate(mols):
            post_adj = mol.GetAdjacencyMatrix()
            try:
                (pre_adj == post_adj).all()
            except:
                print(mol_confs_sdf)
                break
            
            if (pre_adj == post_adj).all():
                valid_mol[mol_id][conf_id] = {}
                xyz = mol.ToXYZ()
                en = mol.GetProp("ConfEnergies")
                valid_mol[mol_id][conf_id]["ff_xyz"] = xyz
                valid_mol[mol_id][conf_id]["ff_energy"] = en
            else:
                failed_jobs[mol_id][conf_id] = 'adjacency matrix'
            
        return failed_jobs, valid_mol
    else:
        return None

input_smiles_path = sys.argv[1]
output_file_name = sys.argv[2]
n_jobs = int(sys.argv[3])

submit_dir = os.getcwd()

# input_smiles_path = "reactants_products_wb97xd_and_xtb_opted_ts_combo_results_hashed_chart_aug11b.csv"
# n_jobs = 8

df = pd.read_csv(input_smiles_path)
mol_ids = list(df.id)
mol_id_to_smi = dict(zip(df.id, df.smiles))

out = Parallel(n_jobs=n_jobs, backend="multiprocessing", verbose=5)(delayed(parser)(mol_id) for mol_id in mol_ids)
out = [x for x in out if x is not None]

with open(os.path.join(submit_dir, f'{output_file_name}.pkl'), 'wb') as outfile:
    pkl.dump(out, outfile)

xyz_FF_opt = {}
for failed_dict, success_dict in out:
    for mol_id in success_dict:
        xyz_FF_opt[mol_id] = {conf_id: success_dict[mol_id][conf_id]["ff_xyz"] for conf_id in success_dict[mol_id]}

with open(f"{output_file_name}_xyz.pkl", "wb") as f:
    pkl.dump(xyz_FF_opt, f)