#!/usr/bin/env python
# coding: utf-8
import os
import sys
import pandas as pd
import pickle as pkl
from tqdm import tqdm
from joblib import Parallel, delayed

from radical_workflow.parser.dft_opt_freq_parser import dft_opt_freq_parser

input_smiles_path = sys.argv[1]
output_file_name = sys.argv[2]
n_jobs = int(sys.argv[3])

# input_smiles_path = "reactants_products_wb97xd_and_xtb_opted_ts_combo_results_hashed_chart_aug11b.csv"
# n_jobs = 8

df = pd.read_csv(input_smiles_path)
mol_ids = df['id'].tolist()
mol_id_to_smi = dict(zip(df['id'], df['smiles']))

out = Parallel(n_jobs=n_jobs, backend="multiprocessing", verbose=5)(delayed(dft_opt_freq_parser)(mol_id, mol_id_to_smi[mol_id]) for mol_id in tqdm(mol_ids))

failed_jobs = dict()
valid_jobs = dict()
for failed_job, valid_job in out:
    failed_jobs.update(failed_job)
    valid_jobs.update(valid_job)

with open(os.path.join(f'{output_file_name}.pkl'), 'wb') as outfile:
    pkl.dump(valid_jobs, outfile, protocol=pkl.HIGHEST_PROTOCOL)

with open(os.path.join(f'{output_file_name}_failed.pkl'), 'wb') as outfile:
    pkl.dump(failed_jobs, outfile, protocol=pkl.HIGHEST_PROTOCOL)

mol_id_to_DFT_opted_xyz_std_ori = {}
for mol_id, valid_job in valid_jobs.items():
    mol_id_to_DFT_opted_xyz_std_ori[mol_id] = valid_job['dft_xyz_std_ori']

with open(os.path.join(f'{output_file_name}_xyz_std_ori.pkl'), 'wb') as outfile:
    pkl.dump(mol_id_to_DFT_opted_xyz_std_ori, outfile, protocol=pkl.HIGHEST_PROTOCOL)

mol_id_to_DFT_opted_xyz = {}
for mol_id, valid_job in valid_jobs.items():
    mol_id_to_DFT_opted_xyz[mol_id] = valid_job['dft_xyz']

with open(os.path.join(f'{output_file_name}_xyz.pkl'), 'wb') as outfile:
    pkl.dump(mol_id_to_DFT_opted_xyz, outfile, protocol=pkl.HIGHEST_PROTOCOL)

print(f"Total number of molecules: {len(mol_ids)}")
print(f"Total number of failed molecules: {len(failed_jobs)}")
print(failed_jobs)

print("Done!")