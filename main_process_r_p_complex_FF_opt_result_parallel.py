import os
import sys
import csv
import tarfile
import pickle as pkl
import pandas as pd
from joblib import Parallel, delayed

# import sys
# sys.path.insert(0, "/home/gridsan/hwpang/RMG_shared/Software/RDMC-main/")

from rdkit import Chem
from rdmc.mol import RDKitMol

def parser(ts_id, submit_dir):
    ids = str(int(int(ts_id.split("id")[1])/1000))
    tar_file_path = os.path.join(submit_dir, "output", "r_p_complex_FF_opt", "outputs", f"outputs_{ids}", f"{ts_id}.tar")
    if os.path.exists(tar_file_path):
        tar = tarfile.open(tar_file_path)
        ts_smi = ts_id_to_smi[ts_id]
        r_smi, p_smi = ts_smi.split(">>")
        ts_xyz = ts_id_to_xyz[ts_id]
        ts_mol = RDKitMol.FromXYZ(ts_xyz, header=False, sanitize=False)
        ts_mol = ts_mol._mol
        ts_mol.SetProp("_Name", ts_id)
        for member in tar:
            if "_r.sdf" in member.name:
                f = tar.extractfile(member)
                mols = Chem.ForwardSDMolSupplier(f, removeHs=False)
                r_mol = [mol for mol in mols][0]
            elif "_p.sdf" in member.name:
                f = tar.extractfile(member)
                mols = Chem.ForwardSDMolSupplier(f, removeHs=False)
                p_mol = [mol for mol in mols][0]    
        tar.close()
        return ts_id, r_smi, p_smi, ts_mol, r_mol, p_mol
    else:
        return None

input_smiles_path = sys.argv[1]
output_file_name = sys.argv[2]
n_jobs = int(sys.argv[3])
ts_id_to_xyz_path = sys.argv[4]

submit_dir = os.getcwd()

# input_smiles_path = "/home/gridsan/groups/RMG/Projects/Hao-Wei-Oscar-Yunsie/production_run/HAbs/inputs/TS_sep1a_all/wb97xd_and_xtb_opted_ts_combo_results_hashed_sep1a_ts_input.csv"
# n_jobs = 8
# output_file_name = "test"
# ts_id_to_xyz_path = "/home/gridsan/groups/RMG/Projects/Hao-Wei-Oscar-Yunsie/production_run/HAbs/inputs/TS_sep1a_all/wb97xd_and_xtb_opted_ts_combo_results_hashed_sep1a_ts_dft_xyz.pkl"

df = pd.read_csv(input_smiles_path)
ts_id_to_smi = dict(zip(df.id, df.rxn_smiles))
ts_ids = list(df.id)
# ts_ids = list(df.id)[:10]
with open(ts_id_to_xyz_path, "rb") as f:
    ts_id_to_xyz = pkl.load(f)

out = Parallel(n_jobs=n_jobs, backend="multiprocessing", verbose=5)(delayed(parser)(ts_id, submit_dir) for ts_id in ts_ids)
out = [x for x in out if x is not None]
csv_file = os.path.join(submit_dir, f'{output_file_name}.csv')
ts_writer = Chem.rdmolfiles.SDWriter(os.path.join(submit_dir, f'{output_file_name}_ts.sdf'))
r_writer = Chem.rdmolfiles.SDWriter(os.path.join(submit_dir, f'{output_file_name}_reactants.sdf'))
p_writer = Chem.rdmolfiles.SDWriter(os.path.join(submit_dir, f'{output_file_name}_products.sdf'))

header = ['id', 'rsmi', 'psmi']

with open(csv_file , 'w') as csvfile:
    # creating a csv writer object
    csvwriter = csv.writer(csvfile)
    # writing the header
    csvwriter.writerow(header)

    for ts_id, r_smi, p_smi, ts_mol, r_mol, p_mol in out:
        csvwriter.writerow([ts_id, r_smi, p_smi])
        r_mol.SetProp("_Name", r_smi)
        p_mol.SetProp("_Name", p_smi)
        ts_writer.write(ts_mol)
        r_writer.write(r_mol)
        p_writer.write(p_mol)

ts_writer.close()
r_writer.close()
p_writer.close()