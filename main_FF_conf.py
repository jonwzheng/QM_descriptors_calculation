from argparse import ArgumentParser
import os
import pandas as pd
from lib.FF_conf_generation import _genConf

parser = ArgumentParser()
parser.add_argument('--input_smiles', type=str, required=True,
                    help='input smiles included in a .csv file')
parser.add_argument('--output_folder', type=str, default='output',
                    help='output folder name')

# conformer searching
parser.add_argument('--FF_conf_folder', type=str, default='FF_conf',
                    help='Folder name for FF searched conformers')
parser.add_argument('--max_n_conf', type=int, default=800,
                    help='maximum number of FF conformers. nc = 3**n_rotatable_bonds, n_conf = nc if nc < max_n_conf else max_n_conf')
parser.add_argument('-max_conf_try', type=int, default=2000,
                    help='maximum attempt for conformer generating, '
                         'this is useful for molecules with many chiral centers.')
parser.add_argument('-rmspre', type=float, required=False, default=0.1,
                        help='rms threshold pre optimization')
parser.add_argument('--rmspost', type=float, required=False, default=0.4,
                    help='rms threshold post FF minimization')
parser.add_argument('--E_cutoff_fraction', type=float, required=False, default=0.5,
                    help='energy window for FF minimization.')
parser.add_argument('--n_lowest_E_confs_to_save', type=int, default=10,
                    help='number of lowest energy conformers to save')

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
parser.add_argument('--scratch_dir', type=str, required=True,
                    help='scratch directory')

args = parser.parse_args()

XTB_PATH = args.XTB_path
G16_PATH = args.G16_path
RDMC_PATH = args.RDMC_path
COSMOTHERM_PATH = args.COSMOtherm_path
COSMO_DATABASE_PATH = args.COSMO_database_path
ORCA_PATH = args.ORCA_path

submit_dir = os.path.abspath(os.getcwd())
output_dir = os.path.join(submit_dir, args.output_folder)
FF_conf_dir = os.path.join(output_dir, args.FF_conf_folder)

df = pd.read_csv(args.input_smiles, index_col=0)
assert len(df['id']) == len(set(df['id'])), "ids must be unique"

# conformer searching
conf_search_FFs = ["GFNFF", "MMFF94s"]

assert XTB_PATH is not None, f"XTB_PATH must be provided to use GFNFF"

mol_ids = list(df["id"])
smiles_list = list(df["smiles"])
mol_id_to_smiles = dict(zip(mol_ids, smiles_list))

os.makedirs(args.scratch_dir, exist_ok=True)

for conf_search_FF in conf_search_FFs:
    for _ in range(5):
        for subinputs_folder in os.listdir(os.path.join(FF_conf_dir, "inputs")):
            ids = subinputs_folder.split("_")[1]
            subinputs_dir = os.path.join(FF_conf_dir, "inputs", subinputs_folder)
            suboutputs_dir = os.path.join(FF_conf_dir, "outputs", f"outputs_{ids}")
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
                        smi = mol_id_to_smiles[mol_id]
                        print(smi)
                        _genConf(smi, mol_id, XTB_PATH, conf_search_FF, args.max_n_conf, args.max_conf_try, args.rmspre, args.E_cutoff_fraction, args.rmspost, args.n_lowest_E_confs_to_save, args.scratch_dir, suboutputs_dir, subinputs_dir)
