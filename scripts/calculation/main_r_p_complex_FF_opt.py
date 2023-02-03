from argparse import ArgumentParser
import os
import pickle as pkl
import pandas as pd
from rdkit import Chem

from lib.reset_r_p_complex import reset_r_p_complex_FF_opt

parser = ArgumentParser()
parser.add_argument('--input_smiles', type=str, required=True,
                    help='input smiles included in a .csv file')
parser.add_argument('--output_folder', type=str, default='output',
                    help='output folder name')
parser.add_argument('--scratch_dir', type=str, required=True,
                    help='scratch dir')
parser.add_argument('--xyz_DFT_opt_dict', type=str, required=True,
                    help='pickled dict mapping from ts_id to xyz')

# reactant complex and product complex semiempirical optimization calculation
parser.add_argument('--r_p_complex_FF_opt_folder', type=str, default='r_p_complex_FF_opt',
                    help='folder for reactant complex and product complex semiempirical optimization')

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
r_p_complex_FF_opt_dir = os.path.join(output_dir, args.r_p_complex_FF_opt_folder)

inputs_dir = os.path.join(r_p_complex_FF_opt_dir, "inputs")
outputs_dir = os.path.join(r_p_complex_FF_opt_dir, "outputs")

df = pd.read_csv(args.input_smiles, index_col=0)
assert len(df['id']) == len(set(df['id'])), "ids must be unique"

with open(args.xyz_DFT_opt_dict, "rb") as f:
    xyz_DFT_opt_dict = pkl.load(f)

assert RDMC_PATH is not None, "RDMC_PATH must be provided for FF opt"

ts_ids = list(df["id"])
rxn_smiles_list = list(df["rxn_smiles"])
ts_id_to_rxn_smi = dict(zip(ts_ids, rxn_smiles_list))
os.makedirs(args.scratch_dir, exist_ok=True)

for _ in range(5):
    for subinputs_folder in os.listdir(inputs_dir):
        ids = subinputs_folder.split("_")[1]
        subinputs_dir = os.path.join(inputs_dir, subinputs_folder)
        suboutputs_dir = os.path.join(outputs_dir, f"outputs_{ids}")
        for input_file in os.listdir(subinputs_dir):
            if ".in" in input_file:
                ts_id = input_file.split(".in")[0]
                print(ts_id)
                try:
                    os.rename(os.path.join(subinputs_dir, input_file), os.path.join(subinputs_dir, f"{ts_id}.tmp"))
                except:
                    continue
                else:
                    rxn_smi = ts_id_to_rxn_smi[ts_id]
                    ts_xyz = xyz_DFT_opt_dict[ts_id]
                    reset_r_p_complex_FF_opt(rxn_smi, ts_xyz, ts_id, subinputs_dir, suboutputs_dir, args.scratch_dir)
