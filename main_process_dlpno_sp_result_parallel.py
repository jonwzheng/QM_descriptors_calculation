import os
import sys
import pickle as pkl
import pandas as pd
from joblib import Parallel, delayed

class OrcaLog(object):
    def __init__(self, path):
        self.path = path

    def check_for_errors(self):
        """
        Checks for common errors in an Orca log file.
        If any are found, this method will raise an error and crash.
        """
        with open(os.path.join(self.path), 'r') as f:
            lines = f.readlines()
            error = None
            for line in reversed(lines):
                # check for common error messages
                if 'ORCA finished by error termination in SCF' in line:
                    error = 'SCF'
                    break
                elif 'ORCA finished by error termination in MDCI' in line:
                    error = 'MDCI'
                    break
                elif 'Error : multiplicity' in line:
                    error = f'The multiplicity and charge combination are wrong.'
                    break
                elif 'ORCA TERMINATED NORMALLY' in line:
                    break
        return error

    def load_energy(self):
        """
        Load the energy in J/ml from an Orca log file. Only the last energy
        in the file is returned. The zero-point energy is *not* included in
        the returned value.
        """
        e_elect = None
        with open(self.path, 'r') as f:
            for line in f:
                if 'FINAL SINGLE POINT ENERGY' in line:  # for all methods in Orca
                    e_elect = float(line.split()[-1])
        if e_elect is None:
            raise LogError('Unable to find energy in Orca output file.')
        return e_elect

def parser(mol_log):

    failed_jobs = dict()
    valid_mol = dict()


    mol_id = os.path.basename(mol_log).split(".log")[0]
    mol_smi = df.loc[df['id'] == mol_id]['smiles'].tolist()[0]

    orca_log = mol_log

    olog = OrcaLog(orca_log)

    error = olog.check_for_errors()

    if error is not None:
        failed_jobs[mol_id] = dict()
        failed_jobs[mol_id]['status'] = False
        failed_jobs[mol_id]['reason'] = error
        return failed_jobs, valid_mol

    valid_mol[mol_id] = dict()
    valid_mol[mol_id]['mol_smi'] = mol_smi
    valid_mol[mol_id]['dlpno_energy'] = olog.load_energy()

    return failed_jobs, valid_mol

input_smiles_path = sys.argv[1]
output_file_name = sys.argv[2]
n_jobs = int(sys.argv[3])

df = pd.read_csv(input_smiles_path)
mol_log_paths = []
submit_dir = os.getcwd()
for suboutput_folder in os.listdir(os.path.join(submit_dir, "output", "DLPNO_sp", "outputs")):
    for mol_log in os.listdir(os.path.join(submit_dir, "output", "DLPNO_sp", "outputs", suboutput_folder)):
        if ".log" in mol_log:
            mol_log_paths.append(os.path.join(submit_dir, "output", "DLPNO_sp", "outputs", suboutput_folder, mol_log))

out = Parallel(n_jobs=n_jobs, backend="multiprocessing", verbose=5)(delayed(parser)(mol_log) for mol_log in mol_log_paths)

with open(os.path.join(submit_dir, f'{output_file_name}.pkl'), 'wb') as outfile:
    pkl.dump(out, outfile)