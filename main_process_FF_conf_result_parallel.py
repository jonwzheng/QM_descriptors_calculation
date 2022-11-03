#!/usr/bin/env python
# coding: utf-8

from argparse import ArgumentParser

import os
import re
import sys
import shutil

import numpy as np
import pandas as pd
import pickle as pkl

from lib.file_parser import load_sdf

from joblib import Parallel, delayed

def parser(mol_confs_sdf):

    failed_jobs = dict()
    valid_mol = dict()

    mol_id = os.path.basename(mol_confs_sdf).split("_confs.sdf")[0]
    mol_smi = df.loc[df['id'] == mol_id]['smiles'].tolist()[0]
    pre_adj = Chem.GetAdjacencyMatrix(Chem.MolFromSmiles(mol_smi))

    opt_mol = load_sdf(mol_confs_sdf)
    conf_ids = list(range(opt_mol.GetNumConformers()))
    for conf_id in conf_ids:
        xyz = Chem.MolToXYZBlock(opt_mol, conf_id)
        en = opt_mol.ConfEnergies
    return failed_jobs, valid_mol

input_smiles_path = sys.argv[1]
output_file_name = sys.argv[2]
n_jobs = int(sys.argv[3])

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