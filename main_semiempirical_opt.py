from argparse import ArgumentParser
import os
import pickle as pkl
import pandas as pd
from rdkit import Chem

from lib.semiempirical_calculation import semiempirical_opt

parser = ArgumentParser()
parser.add_argument('--input_smiles', type=str, required=True,
                    help='input smiles included in a .csv file')
parser.add_argument('--output_folder', type=str, default='output',
                    help='output folder name')
parser.add_argument('--scratch_dir', type=str, required=True,
                    help='scratch dir')
parser.add_argument('--xyz_FF_dict', type=str, required=True,
                    help='pickled dict mapping from mol_id to confs xyz')

# semiempirical optimization calculation
parser.add_argument('--semiempirical_opt_folder', type=str, default='semiempirical_opt',
                    help='folder for semiempirical optimization')
parser.add_argument('--gaussian_semiempirical_opt_theory', type=str, default='#opt=(calcall,maxcycle=128,noeig,nomicro,cartesian)',
                    help='level of theory for the Gaussian semiempirical calculation')
parser.add_argument('--gaussian_semiempirical_opt_n_procs', type=int, default=4,
                    help='number of process for Gaussian semiempirical calculations')
parser.add_argument('--gaussian_semiempirical_opt_job_ram', type=int, default=2000,
                    help='amount of ram (MB) allocated for each Gaussian semiempirical calculation')

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
semiempirical_opt_dir = os.path.join(output_dir, args.semiempirical_opt_folder)

df = pd.read_csv(args.input_smiles, index_col=0)
assert len(df['id']) == len(set(df['id'])), "ids must be unique"

with open(args.xyz_FF_dict, "rb") as f:
    xyz_FF_dict = pkl.load(f)

assert XTB_PATH is not None, "XTB_PATH must be provided for semiempirical opt"
assert G16_PATH is not None, "G16_PATH must be provided for semiempirical opt"
assert RDMC_PATH is not None, "RDMC_PATH must be provided for semiempirical opt"

mol_ids = list(df["id"])
smiles_list = list(df["smiles"])

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

os.makedirs(args.scratch_dir, exist_ok=True)

semiempirical_opt_methods = ["GFN2-XTB", "pm7", "am1"]

for semiempirical_opt_method in semiempirical_opt_methods:
    for _ in range(5):
        for subinputs_folder in os.listdir(os.path.join(semiempirical_opt_dir, "inputs")):
            ids = subinputs_folder.split("_")[1]
            subinputs_dir = os.path.join(semiempirical_opt_dir, "inputs", subinputs_folder)
            suboutputs_dir = os.path.join(semiempirical_opt_dir, "outputs", f"outputs_{ids}")
            for input_file in os.listdir(subinputs_dir):
                if ".in" in input_file:
                    mol_id = input_file.split(".in")[0]
                    print(mol_id)
                    try:
                        os.rename(os.path.join(subinputs_dir, input_file), os.path.join(subinputs_dir, f"{mol_id}.tmp"))
                    except:
                        continue
                    else:
                        ids = str(int(int(mol_id.split("id")[1])/1000))
                        smi = mol_id_to_smi_dict[mol_id]
                        charge = mol_id_to_charge_dict[mol_id]
                        mult = mol_id_to_mult_dict[mol_id]
                        print(smi)
                        tmp_mol_dir = os.path.join(subinputs_dir, mol_id)
                        os.makedirs(tmp_mol_dir, exist_ok=True)
                        semiempirical_opt(mol_id, charge, mult, xyz_FF_dict, XTB_PATH, RDMC_PATH, G16_PATH, args.gaussian_semiempirical_opt_theory, args.gaussian_semiempirical_opt_n_procs, args.gaussian_semiempirical_opt_job_ram, semiempirical_opt_method, args.scratch_dir, tmp_mol_dir, suboutputs_dir, subinputs_dir)
