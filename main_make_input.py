from argparse import ArgumentParser
import os
import shutil
import time
import tarfile
import csv

import pickle as pkl
import pandas as pd
import traceback

import rdkit.Chem as Chem

from lib.wft_calculation import generate_dlpno_sp_input

parser = ArgumentParser()
parser.add_argument('--input_smiles', type=str, required=True,
                    help='input smiles included in a .csv file')
parser.add_argument('--output_folder', type=str, default='output',
                    help='output folder name')
parser.add_argument('--xyz_DFT_opt', type=str, default=None,
                    help='pickle file containing a dictionary to map between the mol_id and DFT-optimized xyz for following calculations',)
parser.add_argument('--DLPNO_sp_folder', type=str, default='DLPNO_sp')
parser.add_argument('--DLPNO_sp_n_procs', type=int, default=22,
                    help='number of process for DLPNO calculations')
parser.add_argument('--DLPNO_sp_job_ram', type=int, default=3900,
                    help='amount of ram (MB) per core allocated for each DLPNO calculation')

args = parser.parse_args()

# input files
with open(args.xyz_DFT_opt, "rb") as f:
    xyz_DFT_opt = pkl.load(f)

df = pd.read_csv(args.input_smiles, index_col=0)

# create id to smile mapping
mol_id_to_smi_dict = dict(zip(df.id, df.smiles))
mol_id_to_charge_dict = dict()
mol_id_to_mult_dict = dict()
for k, v in mol_id_to_smi_dict.items():
    try:
        mol = Chem.MolFromSmiles(v)
    except Exception as e:
        print(f'Cannot translate smi {v} to molecule for species {k}')

    try:
        charge = Chem.GetFormalCharge(mol)
        mol_id_to_charge_dict[k] = charge
    except Exception as e:
        print(f'Cannot determine molecular charge for species {k} with smi {v}')

    num_radical_elec = 0
    for atom in mol.GetAtoms():
        num_radical_elec += atom.GetNumRadicalElectrons()
    mol_id_to_mult_dict[k] =  num_radical_elec + 1

submit_dir = os.path.abspath(os.getcwd())
project_dir = os.path.abspath(os.path.join(args.output_folder))

os.makedirs(os.path.join(project_dir, args.DLPNO_sp_folder), exist_ok=True)
os.makedirs(os.path.join(project_dir, args.DLPNO_sp_folder, "inputs"), exist_ok=True)
os.makedirs(os.path.join(project_dir, args.DLPNO_sp_folder, "outputs"), exist_ok=True)

mol_ids = list(df["id"])
mol_ids = [mol_id for mol_id in mol_ids if mol_id in xyz_DFT_opt]
for mol_id in mol_ids:
    ids = str(int(int(mol_id.split("id")[1])/1000))
    os.makedirs(os.path.join(project_dir, args.DLPNO_sp_folder, "inputs", f"inputs_{ids}"), exist_ok=True)
    os.makedirs(os.path.join(project_dir, args.DLPNO_sp_folder, "outputs", f"outputs_{ids}"), exist_ok=True)
    if not os.path.exists(os.path.join(project_dir, args.DLPNO_sp_folder, "inputs", f"inputs_{ids}", f"{mol_id}.in")):
        charge = mol_id_to_charge_dict[mol_id]
        mult = mol_id_to_mult_dict[mol_id]
        coords = xyz_DFT_opt[mol_id]
        script = generate_dlpno_sp_input(coords, charge, mult, args.DLPNO_sp_job_ram, args.DLPNO_sp_n_procs)

        with open(os.path.join(project_dir, args.DLPNO_sp_folder, "inputs", f"inputs_{ids}", f"{mol_id}.in"), "w+") as f:
            f.write(script)
    else:
        continue