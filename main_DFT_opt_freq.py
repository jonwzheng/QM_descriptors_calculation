from argparse import ArgumentParser
import os
import pickle as pkl
import pandas as pd
from rdkit import Chem

from lib.dft_calculation import dft_scf_opt

parser = ArgumentParser()
parser.add_argument('--input_smiles', type=str, required=True,
                    help='input smiles included in a .csv file')
parser.add_argument('--output_folder', type=str, default='output',
                    help='output folder name')
parser.add_argument('--scratch_dir', type=str, required=True,
                    help='scratch dir')
parser.add_argument('--xyz_semiempirical_opt_dict', type=str, required=True,
                    help='pickled dict mapping from mol_id to semiempirical opted xyz')

# DFT optimization and frequency calculation
parser.add_argument('--DFT_opt_freq_folder', type=str, default='DFT_opt_freq',
                    help='folder for DFT optimization and frequency calculation',)
parser.add_argument('--DFT_opt_freq_theory', type=str, default='#P opt=(calcfc,maxcycle=128,noeig,nomicro,cartesian) freq scf=(xqc) iop(7/33=1) iop(2/9=2000) guess=mix wb97xd/def2svp',
                    help='level of theory for the DFT calculation')
parser.add_argument('--DFT_opt_freq_theory_backup', type=str, default='#P opt=(calcall,maxcycle=64,noeig,nomicro,cartesian) freq scf=(tight, xqc) iop(7/33=1) iop(2/9=2000) guess=mix wb97xd/def2svp',
                    help='level of theory for the DFT calculation if DFT_opt_freq_theory failed')
parser.add_argument('--DFT_opt_freq_n_procs', type=int, default=16,
                    help='number of process for DFT calculations')
parser.add_argument('--DFT_opt_freq_job_ram', type=int, default=62400, #3900*16
                    help='amount of ram (MB) allocated for each DFT calculation')

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
DFT_opt_freq_dir = os.path.join(output_dir, args.DFT_opt_freq_folder)

df = pd.read_csv(args.input_smiles, index_col=0)
assert len(df['id']) == len(set(df['id'])), "ids must be unique"

with open(args.xyz_semiempirical_opt_dict, "rb") as f:
    xyz_semiempirical_opt_dict = pkl.load(f)

assert G16_PATH is not None, "G16_PATH must be provided for DFT opt and freq"

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

DFT_opt_freq_theories = [args.DFT_opt_freq_theory, args.DFT_opt_freq_theory_backup]

for DFT_opt_freq_theory in DFT_opt_freq_theories:
    for _ in range(5):
        for subinputs_folder in os.listdir(os.path.join(DFT_opt_freq_dir, "inputs")):
            ids = subinputs_folder.split("_")[1]
            subinputs_dir = os.path.join(DFT_opt_freq_dir, "inputs", subinputs_folder)
            suboutputs_dir = os.path.join(DFT_opt_freq_dir, "outputs", f"outputs_{ids}")
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
                        dft_scf_opt(mol_id, xyz_semiempirical_opt_dict, G16_PATH, DFT_opt_freq_theory, args.DFT_opt_freq_n_procs, args.DFT_opt_freq_job_ram, charge, mult, args.scratch_dir, suboutputs_dir, subinputs_dir)