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
parser.add_argument('--xyz_DFT_opt_dict', type=str, default=None,
                    help='pickle file containing a dictionary to map between the mol_id and DFT-optimized xyz for following calculations',)

#dlpno sp
parser.add_argument('--DLPNO_sp_folder', type=str, default='DLPNO_sp')
parser.add_argument('--DLPNO_sp_n_procs', type=int, default=22,
                    help='number of process for DLPNO calculations')
parser.add_argument('--DLPNO_sp_job_ram', type=int, default=3900,
                    help='amount of ram (MB) per core allocated for each DLPNO calculation')

# specify paths
parser.add_argument('--XTB_path', type=str, required=False, default=None,
                    help='path to installed XTB')
parser.add_argument('--G16_path', type=str, required=False, default=None,
                    help='path to installed Gaussian 16')
parser.add_argument('--RDMC_path', type=str, required=False, default=None,
                    help='path to RDMC to use xtb-gaussian script for xtb optimization calculation.')
parser.add_argument('--COSMOtherm_path', type=str, required=False, default=None,
                    help='path to COSMOthermo')
parser.add_argument('--COSMO_database_path', type=str, required=False, default=None,
                    help='path to COSMO_database')
parser.add_argument('--ORCA_path', type=str, required=False, default=None,
                    help='path to ORCA')

args = parser.parse_args()

XTB_PATH = args.XTB_path
G16_PATH = args.G16_path
RDMC_PATH = args.RDMC_path
COSMOTHERM_PATH = args.COSMOtherm_path
COSMO_DATABASE_PATH = args.COSMO_database_path
ORCA_PATH = args.ORCA_path

submit_dir = os.path.abspath(os.getcwd())
output_dir = os.path.join(submit_dir, args.output_folder)
DLPNO_sp_dir = os.path.join(output_dir, args.DLPNO_sp_folder)
os.makedirs(DLPNO_sp_dir, exist_ok=True)

df = pd.read_csv(args.input_smiles, index_col=0)
assert len(df['id']) == len(set(df['id'])), "ids must be unique"

# input files
with open(args.xyz_DFT_opt_dict, "rb") as f:
    xyz_DFT_opt_dict = pkl.load(f)

assert ORCA_PATH is not None, "ORCA_PATH must be provided for dlpno sp calc"

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

mol_ids = list(df["id"])
smiles_list = list(df["smiles"])
inputs_dir = os.path.join(DLPNO_sp_dir, "inputs")
os.makedirs(inputs_dir, exist_ok=True)
outputs_dir = os.path.join(DLPNO_sp_dir, "outputs")
os.makedirs(outputs_dir, exist_ok=True)

for mol_id, smi in zip(mol_ids, smiles_list):
    if mol_id in xyz_DFT_opt_dict:
        ids = str(int(int(mol_id.split("id")[1])/1000))
        subinputs_dir = os.path.join(inputs_dir, f"inputs_{ids}")
        suboutputs_dir = os.path.join(outputs_dir, f"outputs_{ids}")
        os.makedirs(suboutputs_dir, exist_ok=True)
        mol_id_path = os.path.join(subinputs_dir, f"{mol_id}.in")
        if not os.path.exists(os.path.join(suboutputs_dir, f"{mol_id}.log")) and not os.path.exists(mol_id_path):
            os.makedirs(subinputs_dir, exist_ok=True)
            charge = mol_id_to_charge_dict[mol_id]
            mult = mol_id_to_mult_dict[mol_id]
            coords = xyz_DFT_opt_dict[mol_id].strip()
            script = generate_dlpno_sp_input(coords, charge, mult, args.DLPNO_sp_job_ram, args.DLPNO_sp_n_procs)

            with open(mol_id_path, "w+") as f:
                f.write(script)