import os
import sys

sys.path.insert(0, "/home/gridsan/hwpang/RMG_shared/Software/QM_descriptors_calculation-radical_workflow/")
import csv
import tarfile
import pickle as pkl
import pandas as pd
from joblib import Parallel, delayed
from lib.utils import REPLACE_LETTER

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
            temp = line.split(b'T=')[1].split(b'K')[0].strip()  # temp in K

        # get the result values
        if b"Nr Compound" in line:
            line = f.readline()
            solvent_name = line.split()[1]
            line = f.readline()
            solute_name = line.split()[1]
            result_values = line.split()[2:6]  # H (in bar), ln(gamma), pv (vapor pressure in bar), Gsolv (kcal/mol)
            # save the result as one list
            each_data_list.append(
                [solvent_name, solute_name, temp] + result_values + [None])
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
        temp = each_data_list[z][2]
        dGsolv = each_data_list[z][6]
        dGsolv_temp_dict[temp] = dGsolv
        if temp == b'298.15':
            ind_298 = z
    print(dGsolv_temp_dict)
    dGsolv_298 = float(dGsolv_temp_dict[b'298.15'])
    dSsolv_298 = - (float(dGsolv_temp_dict[b'299.15']) - float(dGsolv_temp_dict[b'297.15'])) / (299.15 - 297.15)
    dHsolv_298 = dGsolv_298 + 298.15 * dSsolv_298
    each_data_list[ind_298][7] = '%.8f' % dHsolv_298
    return each_data_list

def parser(tar_file_path):
    each_data_lists = []
    tar = tarfile.open(tar_file_path)
    for member in tar:
        if ".tab" in member.name:
            f = tar.extractfile(member)
            each_data_list = read_cosmo_tab_result_from_tar(f)
            each_data_list = get_dHsolv_value(each_data_list)
    tar.close()
    return each_data_lists

input_smiles_path = sys.argv[1]
output_file_name = sys.argv[2]
n_jobs = int(sys.argv[3])

df = pd.read_csv(input_smiles_path)
tar_file_paths = []
submit_dir = os.getcwd()
for suboutput_folder in os.listdir(os.path.join(submit_dir, "output", "COSMO_calc", "outputs")):
    for tar_file in os.listdir(os.path.join(submit_dir, "output", "COSMO_calc", "outputs", suboutput_folder)):
        if ".tar" in tar_file:
            tar_file_paths.append(os.path.join(submit_dir, "output", "COSMO_calc", "outputs", suboutput_folder, tar_file))

out = Parallel(n_jobs=n_jobs, backend="multiprocessing", verbose=5)(delayed(parser)(tar_file) for tar_file in tar_file_paths)

csv_file = os.path.join(submit_dir, f'{output_file_name}.csv')

header = ['solvent_name', 'solute_name', 'temp (K)',
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
    pkl.dump(df_result, outfile)

os.remove(csv_file)