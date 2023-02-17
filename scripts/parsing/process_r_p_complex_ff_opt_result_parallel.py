import os
import sys
import csv
import tarfile
import pickle as pkl
import pandas as pd
from joblib import Parallel, delayed

from rdkit import Chem
from rdmc.mol import RDKitMol

def parser(ts_id):
    ids = str(int(int(ts_id.split("id")[1])/1000))
    ts_smi = ts_id_to_smi[ts_id]
    r_smi, p_smi = ts_smi.split(">>")

    pre_r_mol = RDKitMol.FromSmiles(r_smi, removeHs=False, sanitize=True)
    pre_p_mol = RDKitMol.FromSmiles(p_smi, removeHs=False, sanitize=True)

    pre_r_adj = pre_r_mol.GetAdjacencyMatrix()
    pre_p_adj = pre_p_mol.GetAdjacencyMatrix()

    ts_xyz = ts_id_to_xyz[ts_id]
    ts_mol = RDKitMol.FromXYZ(ts_xyz, header=False, sanitize=False)
    ts_mol._mol.SetProp("_Name", ts_id)
    
    tar_file_path = os.path.join("output", "r_p_complex_ff_opt", "outputs", f"outputs_{ids}", f"{ts_id}.tar")
    sdf_file_path = os.path.join("output", "r_p_complex_ff_opt", "outputs", f"outputs_{ids}", f"{ts_id}_r_p.sdf")
    if os.path.exists(tar_file_path):
        tar = tarfile.open(tar_file_path)
        for member in tar:
            if "_r.sdf" in member.name:
                f = tar.extractfile(member)
                mols = Chem.ForwardSDMolSupplier(f, removeHs=False, sanitize=True)
                r_mol = [mol for mol in mols][0]
                r_mol = RDKitMol.FromMol(r_mol)

            elif "_p.sdf" in member.name:
                f = tar.extractfile(member)
                mols = Chem.ForwardSDMolSupplier(f, removeHs=False, sanitize=True)
                p_mol = [mol for mol in mols][0]  
                p_mol = RDKitMol.FromMol(p_mol)
        tar.close()
    elif os.path.exists(sdf_file_path):
        r_mol, p_mol = Chem.ForwardSDMolSupplier(sdf_file_path, removeHs=False, sanitize=True)
        r_mol = RDKitMol.FromMol(r_mol)
        p_mol = RDKitMol.FromMol(p_mol)
    else:
        return None

    if any(atom.GetFormalCharge() != 0 for atom in r_mol.GetAtoms()):
        return None
    if any(atom.GetFormalCharge() != 0 for atom in p_mol.GetAtoms()):
        return None

    r_adj = r_mol.GetAdjacencyMatrix()
    if not (r_adj == pre_r_adj).all():
        return None
    p_adj = p_mol.GetAdjacencyMatrix()
    if not (p_adj == pre_p_adj).all():
        return None
    
    return ts_id, r_smi, p_smi, ts_mol._mol, r_mol._mol, p_mol._mol

input_smiles_path = sys.argv[1]
output_file_name = sys.argv[2]
n_jobs = int(sys.argv[3])
try:
    ts_id_to_xyz_path = sys.argv[4]
except IndexError:
    ts_id_to_xyz_path = None

df = pd.read_csv(input_smiles_path)
if "rxn_smi" in df.columns:
    rxn_smis = df.rxn_smi
elif "rxn_smiles" in df.columns:
    rxn_smis = df.rxn_smiles
else:
    raise ValueError("No reaction smiles provided")
ts_id_to_smi = dict(zip(df.id, rxn_smis))
ts_ids = list(df.id)
if ts_id_to_xyz_path is not None:
    with open(ts_id_to_xyz_path, "rb") as f:
        ts_id_to_xyz = pkl.load(f)
elif "dft_xyz" in df.columns:
    ts_id_to_xyz = dict(zip(df.id, df.dft_xyz))
else:
    raise ValueError("No xyz file provided")

out = Parallel(n_jobs=n_jobs, backend="multiprocessing", verbose=5)(delayed(parser)(ts_id) for ts_id in ts_ids)
out = [x for x in out if x is not None]

print(f"Number of TS: {len(ts_ids)}")
print(f"Number of valid TS, r complex, and p complex: {len(out)}")

csv_file = os.path.join(f'{output_file_name}.csv')
ts_writer = Chem.rdmolfiles.SDWriter(os.path.join(f'{output_file_name}_ts.sdf'))
r_writer = Chem.rdmolfiles.SDWriter(os.path.join(f'{output_file_name}_reactants.sdf'))
p_writer = Chem.rdmolfiles.SDWriter(os.path.join(f'{output_file_name}_products.sdf'))

header = ['id', 'rsmi', 'psmi']

with open(csv_file , 'w') as csvfile:
    # creating a csv writer object
    csvwriter = csv.writer(csvfile)
    # writing the header
    csvwriter.writerow(header)

    for ts_id, r_smi, p_smi, ts_mol, r_mol, p_mol in out:
        csvwriter.writerow([ts_id, r_smi, p_smi])
        ts_mol.SetProp("_Name", ts_id)
        r_mol.SetProp("_Name", r_smi)
        p_mol.SetProp("_Name", p_smi)
        ts_writer.write(ts_mol)
        r_writer.write(r_mol)
        p_writer.write(p_mol)

ts_writer.close()
r_writer.close()
p_writer.close()