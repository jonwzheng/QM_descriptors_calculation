from argparse import ArgumentParser
import os
import pandas as pd
import shutil

parser = ArgumentParser()
parser.add_argument('--input_smiles', type=str, required=True,
                    help='input smiles included in a .csv file')
parser.add_argument('--output_folder', type=str, default='output',
                    help='output folder name')
parser.add_argument('--job_id', type=int, default=0,
                    help='slurm job id')
parser.add_argument('--task_id', type=int, default=0,
                    help='task id for job arrays or LLsub')

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
parser.add_argument('--FF_threads', type=int, required=False, default=4,
                    help='number of process for the FF conformer searching')
parser.add_argument('--timeout', required=False, default=7200,
                    help='time window for each FF conformer searching sub process')
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

args = parser.parse_args()

XTB_PATH = args.XTB_path
G16_PATH = args.G16_path
RDMC_PATH = args.RDMC_path
COSMOTHERM_PATH = args.COSMOtherm_path
COSMO_DATABASE_PATH = args.COSMO_database_path
ORCA_PATH = args.ORCA_path

submit_dir = os.path.abspath(os.getcwd())
output_dir = os.path.join(submit_dir, args.output_folder)
os.makedirs(output_dir, exist_ok=True)
FF_conf_dir = os.path.join(output_dir, args.FF_conf_folder)
os.makedirs(FF_conf_dir, exist_ok=True)

df = pd.read_csv(args.input_smiles, index_col=0)
assert len(df['id']) == len(set(df['id'])), "ids must be unique"

assert XTB_PATH is not None, "XTB_PATH must be provided to use GFNFF"

mol_ids = list(df["id"])
smiles_list = list(df["smiles"])
inputs_dir = os.path.join(FF_conf_dir, "inputs")
shutil.rmtree(inputs_dir)
os.makedirs(inputs_dir)
outputs_dir = os.path.join(FF_conf_dir, "outputs")
os.makedirs(outputs_dir, exist_ok=True)

print("Making FF conformer input files...")

for mol_id, smi in zip(mol_ids, smiles_list):
    ids = str(int(int(mol_id.split("id")[1])/1000))
    subinputs_dir = os.path.join(inputs_dir, f"inputs_{ids}")
    suboutputs_dir = os.path.join(outputs_dir, f"outputs_{ids}")
    os.makedirs(suboutputs_dir, exist_ok=True)
    mol_id_path = os.path.join(subinputs_dir, f"{mol_id}.in")
    if not os.path.exists(os.path.join(suboutputs_dir, f"{mol_id}_confs.sdf")) and not os.path.exists(mol_id_path):
        os.makedirs(subinputs_dir, exist_ok=True)
        with open(mol_id_path, "w") as f:
            f.write(mol_id)
        print(mol_id)
    else:
        continue
