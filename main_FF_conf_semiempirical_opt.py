from argparse import ArgumentParser
import os
import pandas as pd
from rdkit import Chem
from rdmc.mol import RDKitMol

from lib.FF_conf_generation import _genConf
from lib.semiempirical_calculation import semiempirical_opt

parser = ArgumentParser()
parser.add_argument('--input_smiles', type=str, required=True,
                    help='input smiles included in a .csv file')
parser.add_argument('--output_folder', type=str, default='output',
                    help='output folder name')
parser.add_argument('--task_id', type=int, default=0,
                    help='task id for the calculation',)
parser.add_argument('--num_tasks', type=int, default=1,
                    help='number of tasks for the calculation',)

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

df = pd.read_csv(args.input_smiles, index_col=0)
assert len(df['id']) == len(set(df['id'])), "ids must be unique"

# conformer searching
conf_search_FFs = ["GFNFF", "MMFF94s"]

assert XTB_PATH is not None, f"XTB_PATH must be provided to use GFNFF"

# create id to smile mapping
mol_ids = df['id'].tolist()
smiles_list = df['smiles'].tolist()
mol_id_to_smi_dict = dict(zip(mol_ids, smiles_list))
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

print("Making input files for conformer searching")

FF_conf_dir = os.path.join(output_dir, args.FF_conf_folder)
inputs_dir = os.path.join(FF_conf_dir, "inputs")
outputs_dir = os.path.join(FF_conf_dir, "outputs")

mol_ids_smis = list(zip(mol_ids, smiles_list))
for mol_id, smi in mol_ids_smis[args.task_id:len(mol_ids_smis):args.num_tasks]:
    ids = str(int(int(mol_id.split("id")[1])/1000))
    subinputs_dir = os.path.join(inputs_dir, f"inputs_{ids}")
    suboutputs_dir = os.path.join(outputs_dir, f"outputs_{ids}")
    os.makedirs(suboutputs_dir, exist_ok=True)
    mol_id_path = os.path.join(subinputs_dir, f"{mol_id}.in")
    if not os.path.exists(os.path.join(suboutputs_dir, f"{mol_id}_confs.sdf")) and not os.path.exists(mol_id_path) and not os.path.join(subinputs_dir, f"{mol_id}.tmp"):
        os.makedirs(subinputs_dir, exist_ok=True)
        with open(mol_id_path, "w") as f:
            f.write(mol_id)
        print(mol_id)

print("Conformer searching with force field...")

FF_conf_mol_ids = []

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
                        smi = mol_id_to_smi_dict[mol_id]
                        print(smi)
                        _genConf(smi, mol_id, XTB_PATH, conf_search_FF, args.max_n_conf, args.max_conf_try, args.rmspre, args.E_cutoff_fraction, args.rmspost, args.n_lowest_E_confs_to_save, args.scratch_dir, suboutputs_dir, subinputs_dir)
                        FF_conf_mol_ids.append(mol_id)


print("Conformer searching with force field done.")

print("Making input files for semiempirical optimization")

semiempirical_opt_dir = os.path.join(output_dir, args.semiempirical_opt_folder)

for mol_id in FF_conf_mol_ids:
    ids = str(int(int(mol_id.split("id")[1])/1000))
    subinputs_dir = os.path.join(semiempirical_opt_dir, "inputs", f"inputs_{ids}")
    os.makedirs(subinputs_dir, exist_ok=True)
    suboutputs_dir = os.path.join(semiempirical_opt_dir, "outputs", f"outputs_{ids}")
    os.makedirs(suboutputs_dir, exist_ok=True)
    if not os.path.exists(os.path.join(subinputs_dir, f"{mol_id}.in")) and not os.path.exists(os.path.join(subinputs_dir, f"{mol_id}.tmp")) and not os.path.exists(os.path.join(suboutputs_dir, f"{mol_id}.tar")):
        with open(os.path.join(subinputs_dir, f"{mol_id}.in"), "w") as f:
            f.write(mol_id)

print("Optimizing conformers with semiempirical method...")

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

                    mol_confs_sdf = os.path.join(FF_conf_dir, "outputs", f"outputs_{ids}", f"{mol_id}_confs.sdf")
                    mols = RDKitMol.FromFile(mol_confs_sdf)
                    xyz_FF_dict = {}
                    xyz_FF_dict[mol_id] = {}
                    for conf_id, mol in enumerate(mols):
                        xyz_FF_dict[mol_id][conf_id] = mol.ToXYZ()
                    
                    try:
                        semiempirical_opt(mol_id, charge, mult, xyz_FF_dict, XTB_PATH, RDMC_PATH, G16_PATH, args.gaussian_semiempirical_opt_theory, args.gaussian_semiempirical_opt_n_procs, args.gaussian_semiempirical_opt_job_ram, args.scratch_dir, tmp_mol_dir, suboutputs_dir, subinputs_dir)
                    except FileNotFoundError as e:
                        print(e)
                        print("Continuing...")
                        continue