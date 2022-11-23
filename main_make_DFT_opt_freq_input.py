from argparse import ArgumentParser
import os
import pickle as pkl
import pandas as pd
import shutil

parser = ArgumentParser()
parser.add_argument('--input_smiles', type=str, required=True,
                    help='input smiles included in a .csv file')
parser.add_argument('--output_folder', type=str, default='output',
                    help='output folder name')
parser.add_argument('--xyz_semiempirical_opt_dict', type=str, required=True,
                    help='pickled dict mapping from mol_id to semiempirical opted xyz')

# DFT optimization and frequency calculation
parser.add_argument('--DFT_opt_freq_folder', type=str, default='DFT_opt_freq',
                    help='folder for DFT optimization and frequency calculation',)
parser.add_argument('--DFT_opt_freq_theory', type=str, default='#P opt=(calcfc,maxcycle=128,noeig,nomicro,cartesian) freq scf=(xqc) iop(7/33=1) iop(2/9=2000) guess=mix wb97xd/def2svp',
                    help='level of theory for the DFT calculation')
parser.add_argument('--DFT_opt_freq_theory_backup', type=str, default='#P opt=(calcall,maxcycle=64,noeig,nomicro,cartesian) freq scf=(tight, xqc) iop(7/33=1) iop(2/9=2000) guess=mix wb97xd/def2svp',
                    help='level of theory for the DFT calculation if DFT_opt_freq_theory failed')
parser.add_argument('--DFT_opt_freq_n_procs', type=int, default=4,
                    help='number of process for DFT calculations')
parser.add_argument('--DFT_opt_freq_job_ram', type=int, default=16000,
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
os.makedirs(DFT_opt_freq_dir, exist_ok=True)

df = pd.read_csv(args.input_smiles, index_col=0)
assert len(df['id']) == len(set(df['id'])), "ids must be unique"

with open(args.xyz_semiempirical_opt_dict, "rb") as f:
    xyz_semiempirical_opt_dict = pkl.load(f)

assert G16_PATH is not None, "G16_PATH must be provided for DFT opt and freq"

mol_ids = list(df["id"])
smiles_list = list(df["smiles"])
inputs_dir = os.path.join(DFT_opt_freq_dir, "inputs")
os.makedirs(inputs_dir, exist_ok=True)
outputs_dir = os.path.join(DFT_opt_freq_dir, "outputs")
os.makedirs(outputs_dir, exist_ok=True)

for mol_id, smi in zip(mol_ids, smiles_list):
    if mol_id in xyz_semiempirical_opt_dict:
        ids = str(int(int(mol_id.split("id")[1])/1000))
        subinputs_dir = os.path.join(inputs_dir, f"inputs_{ids}")
        os.makedirs(subinputs_dir, exist_ok=True)
        suboutputs_dir = os.path.join(outputs_dir, f"outputs_{ids}")
        os.makedirs(suboutputs_dir, exist_ok=True)
        try:
            os.remove(os.path.join(subinputs_dir, f"{mol_id}.tmp"))
        except:
            pass
        mol_id_path = os.path.join(subinputs_dir, f"{mol_id}.in")
        if not os.path.exists(os.path.join(suboutputs_dir, f"{mol_id}.log")) and not os.path.exists(mol_id_path):
            with open(mol_id_path, "w") as f:
                f.write(mol_id)
        else:
            continue
