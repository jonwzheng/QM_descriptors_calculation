import os
import sys
import csv
import tarfile
import pickle as pkl
import pandas as pd
from joblib import Parallel, delayed
from tqdm import tqdm
from radical_workflow.calculation.utils import REPLACE_LETTER

def read_cosmo_tab_result_from_tar(f):
    """
    Modified from Yunsie's code
    """
    each_data_list = []
    # initialize everything
    solvent_name, solute_name, temp = None, None, None
    result_values = None
    line = f.readline()
    while line:
        # get the temperature and mole fraction
        if b"Settings  job" in line:
            temp = line.split(b'T=')[1].split(b'K')[0].strip().decode('utf-8')  # temp in K

        # get the result values
        if b"Nr Compound" in line:
            line = f.readline()
            solvent_name = line.split()[1].decode('utf-8')
            solvent_smiles = solvent_name_to_smi[solvent_name]
            line = f.readline()
            solute_name = line.split()[1].decode('utf-8')
            solute_smiles = mol_id_to_smi[solute_name]
            result_values = line.split()[2:6]  # H (in bar), ln(gamma), pv (vapor pressure in bar), Gsolv (kcal/mol)
            result_values = [result_value.decode('utf-8') for result_value in result_values]
            # save the result as one list
            each_data_list.append(
                [solvent_name, solvent_smiles, solute_name, solute_smiles, temp] + result_values + [None])
            # initialize everything
            solvent_name, solute_name, temp = None, None, None
            result_values = None
        line = f.readline()
    return each_data_list

def get_dHsolv_value(each_data_list):
    # compute solvation enthalpy
    dGsolv_temp_dict = {}
    ind_298 = None
    for z in range(len(each_data_list)):
        temp = each_data_list[z][4]
        dGsolv = each_data_list[z][8]
        dGsolv_temp_dict[temp] = dGsolv
        if temp == '298.15':
            ind_298 = z
    dGsolv_298 = float(dGsolv_temp_dict['298.15'])
    dSsolv_298 = - (float(dGsolv_temp_dict['299.15']) - float(dGsolv_temp_dict['297.15'])) / (299.15 - 297.15)
    dHsolv_298 = dGsolv_298 + 298.15 * dSsolv_298
    each_data_list[ind_298][9] = '%.8f' % dHsolv_298
    return each_data_list

def parser(mol_id):
    ids = str(int(int(mol_id.split("id")[1])/1000))
    tar_file_path = os.path.join(submit_dir, "output", "COSMO_calc", "outputs", f"outputs_{ids}", f"{mol_id}.tar")
    if os.path.isfile(tar_file_path):
        each_data_lists = []
        tar = tarfile.open(tar_file_path)
        for member in tar:
            if ".tab" in member.name:
                f = tar.extractfile(member)
                each_data_list = read_cosmo_tab_result_from_tar(f)
                each_data_list = get_dHsolv_value(each_data_list)
                each_data_lists.append(each_data_list)
        tar.close()
        return each_data_lists
    else:
        return None

input_smiles_path = sys.argv[1]
output_file_name = sys.argv[2]
n_jobs = int(sys.argv[3])
solvent_path = sys.argv[4]

submit_dir = os.getcwd()

# input_smiles_path = "reactants_products_wb97xd_and_xtb_opted_ts_combo_results_hashed_chart_aug11b.csv"
# n_jobs = 8
# output_file_name = "test"

df = pd.read_csv(input_smiles_path)
mol_id_to_smi = dict(zip(df.id, df.smiles))
mol_ids = list(df.id)

df_solvent = pd.read_csv(solvent_path)
solvent_name_to_smi = dict(zip(df_solvent.cosmo_name, df_solvent.smiles))

out = Parallel(n_jobs=n_jobs, backend="multiprocessing", verbose=5)(delayed(parser)(mol_id) for mol_id in tqdm(mol_ids))
failed_mol_ids = [mol_ids[i] for i in range(len(mol_ids)) if out[i] is None]
out = [x for x in out if x is not None]

csv_file = os.path.join(submit_dir, f'{output_file_name}.csv')

header = ['solvent_name', 'solvent_smiles', 'solute_name', 'solute_smiles', 'temp (K)',
        'H (bar)', 'ln(gamma)', 'Pvap (bar)', 'Gsolv (kcal/mol)', 'Hsolv (kcal/mol)']

with open(csv_file , 'w') as csvfile:
    # creating a csv writer object
    csvwriter = csv.writer(csvfile)
    # writing the header
    csvwriter.writerow(header)

    for each_data_lists in out:
        for each_data_list in each_data_lists:
            csvwriter.writerows(each_data_list)

df_result = pd.read_csv(csv_file)

with open(os.path.join(submit_dir, f'{output_file_name}.pkl'), 'wb') as outfile:
    pkl.dump(df_result, outfile, protocol=pkl.HIGHEST_PROTOCOL)

with open(os.path.join(submit_dir, f'{output_file_name}_failed.pkl'), 'wb') as outfile:
    pkl.dump(failed_mol_ids, outfile, protocol=pkl.HIGHEST_PROTOCOL)

print(failed_mol_ids)

os.remove(csv_file)