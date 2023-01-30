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
    failed_job = dict()
    valid_job = dict()

    if os.path.isfile(mol_confs_sdf):

        failed_job[mol_id] = dict()
        valid_job[mol_id] = dict()

        mol_smi = mol_id_to_smi[mol_id]
        pre_adj = RDKitMol.FromSmiles(mol_smi).GetAdjacencyMatrix()
        
        mols = RDKitMol.FromFile(mol_confs_sdf)
        for conf_id, mol in enumerate(mols):
            post_adj = mol.GetAdjacencyMatrix()
            try:
                (pre_adj == post_adj).all()
            except:
                print(mol_confs_sdf)
                break
            
            if (pre_adj == post_adj).all():
                valid_job[mol_id][conf_id] = {}
                xyz = mol.ToXYZ()
                en = mol.GetProp("ConfEnergies")
                valid_job[mol_id][conf_id]["ff_xyz"] = xyz
                valid_job[mol_id][conf_id]["ff_energy"] = en
            else:
                failed_job[mol_id][conf_id] = 'adjacency matrix'
        
        if not valid_job[mol_id]:
            del valid_job[mol_id]
            failed_job[mol_id] = 'all confs failed'
        
        if not failed_job[mol_id]:
            del failed_job[mol_id]
    else:
        failed_job[mol_id] = 'sdf file not found'

    return failed_job, valid_job

input_smiles_path = sys.argv[1]
output_file_name = sys.argv[2]
n_jobs = int(sys.argv[3])

submit_dir = os.getcwd()

##
# input_smiles_path = "inputs/reactants_products_aug11b_inputs.csv"
# output_file_name = "reactants_products_aug11b"
# n_jobs = 1

df = pd.read_csv(input_smiles_path)
mol_ids = list(df.id)
mol_id_to_smi = dict(zip(df.id, df.smiles))

##
# mol_ids = mol_ids[:1000]

out = Parallel(n_jobs=n_jobs, backend="multiprocessing", verbose=5)(delayed(parser)(mol_id) for mol_id in mol_ids)

failed_jobs = dict()
valid_jobs = dict()
for failed_job, valid_job in out:
    failed_jobs.update(failed_job)
    valid_jobs.update(valid_job)

with open(os.path.join(submit_dir, f'{output_file_name}.pkl'), 'wb') as outfile:
    pkl.dump(valid_jobs, outfile)

with open(os.path.join(submit_dir, f'{output_file_name}_failed.pkl'), 'wb') as outfile:
    pkl.dump(failed_jobs, outfile)

xyz_FF_opt = {}
for mol_id in valid_jobs:
    xyz_FF_opt[mol_id] = {conf_id: valid_jobs[mol_id][conf_id]["ff_xyz"] for conf_id in valid_jobs[mol_id]}

with open(f"{output_file_name}_xyz.pkl", "wb") as f:
    pkl.dump(xyz_FF_opt, f)