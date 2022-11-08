import sys
import pickle as pkl
import pandas as pd
from joblib import Parallel, delayed

n_jobs = int(sys.argv[1])

projects = ["aug11b", "sep1a"]

dfs_dict = {}
with open("../inputs/aug11b/reactants_products_wb97xd_and_xtb_opted_ts_combo_results_hashed_chart_aug11b.pkl", "rb") as f:
    dfs_dict["aug11b"] = pkl.load(f)
with open("../inputs/sep1a/reactants_products_wb97xd_and_xtb_opted_ts_combo_results_hashed_lookup_table_sep1a_filtered.pkl", "rb") as f:
    dfs_dict["sep1a"] = pkl.load(f)

def add_project():
    for df_name, df in dfs_dict.items():
        df["project"] = [df_name for i in df.index]
        
add_project()

def create_column_conf(success_dict):
    global prop_names
    global df
    for mol_id in success_dict:
        index = df.index[df["id"]==mol_id][0]
        props_dict = {}
        for prop_name in prop_names:
            props_dict[prop_name] = [success_dict[mol_id][conf_id][prop_name] if conf_id in success_dict[mol_id] else None for conf_id in range(10)]
        
        return index, props_dict

def fill_column_conf(results_dict):
    global prop_names
    global df
    for project, df in dfs_dict.items():
        
        success_dicts = [success_dict for failed_dict, success_dict in results_dict[project] if success_dict]
        
        prop_names = []
        for success_dict in success_dicts:
            for mol_id in success_dict:
                for conf_id in range(10):
                    if conf_id in success_dict[mol_id]:
                        for prop_name in success_dict[mol_id][conf_id]:
                            prop_names.append(prop_name)
                        break
                    break
                break
            break
            
        columns_dict = {}
        for prop_name in prop_names:
            columns_dict.update({f"{prop_name}_conf_{conf_id}": [None for _ in df.index] for conf_id in range(10)})
            
        out = Parallel(n_jobs=n_jobs, backend="multiprocessing", verbose=5)(delayed(create_column_conf)(success_dict) for success_dict in success_dicts)

        for index, props_dict in out:
            for prop_name, props in props_dict.items():
                for conf_id, prop in enumerate(props):
                    columns_dict[f"{prop_name}_conf_{conf_id}"][index] = prop

        for column_name, column in columns_dict.items():
            df[column_name] = column

ff_results_dict = {}
with open("./aug11b/reactants_products_aug11b_ff_opted_results.pkl", "rb") as f:
    ff_results_dict["aug11b"] = pkl.load(f)
with open("./sep1a/reactants_products_sep1a_filtered_ff_opted_results.pkl", "rb") as f:
    ff_results_dict["sep1a"] = pkl.load(f)

fill_column_conf(ff_results_dict)

semi_results_dict = {}
with open("./aug11b/reactants_products_aug11b_semiempirical_opted_results.pkl", "rb") as f:
    semi_results_dict["aug11b"] = pkl.load(f)
with open("./sep1a/reactants_products_sep1a_filtered_semiempirical_opted_results.pkl", "rb") as f:
    semi_results_dict["sep1a"] = pkl.load(f)


fill_column_conf(semi_results_dict)

def create_column(success_dict):
    global prop_names
    global df
    for mol_id in success_dict:
        index = df.index[df["id"]==mol_id][0]
        props_dict = {}
        for prop_name in prop_names:
            props_dict[prop_name] = success_dict[mol_id][prop_name]
        
        return index, props_dict

def fill_column(results_dict):
    global prop_names
    global df
    for project, df in dfs_dict.items():
        
        success_dicts = [success_dict for failed_dict, success_dict in results_dict[project] if success_dict]
        
        prop_names = []
        for success_dict in success_dicts:
            for mol_id in success_dict:
                for prop_name in success_dict[mol_id]:
                    prop_names.append(prop_name)
                break
            break
            
        columns_dict = {}
        for prop_name in prop_names:
            columns_dict.update({f"{prop_name}": [None for _ in df.index]})
            
        out = Parallel(n_jobs=n_jobs, backend="multiprocessing", verbose=5)(delayed(create_column)(success_dict) for success_dict in success_dicts)

        for index, props_dict in out:
            for prop_name, prop in props_dict.items():
                columns_dict[f"{prop_name}"][index] = prop

        for column_name, column in columns_dict.items():
            df[column_name] = column

dft_results_dict = {}
with open("./aug11b/reactants_products_aug11b_dft_opted_results.pkl", "rb") as f:
    dft_results_dict["aug11b"] = pkl.load(f)
with open("./sep1a/reactants_products_sep1a_filtered_dft_opted_results.pkl", "rb") as f:
    dft_results_dict["sep1a"] = pkl.load(f)

fill_column(dft_results_dict)

dlpno_results_dict = {}
with open("./aug11b/reactants_products_aug11b_dlpno_sp_results.pkl", "rb") as f:
    dlpno_results_dict["aug11b"] = pkl.load(f)
with open("./sep1a/reactants_products_sep1a_filtered_dlpno_sp_results.pkl", "rb") as f:
    dlpno_results_dict["sep1a"] = pkl.load(f)

fill_column(dlpno_results_dict)

cosmo_results_dict = {}
with open("./aug11b/reactants_products_aug11b_cosmo_results.pkl", "rb") as f:
    cosmo_results_dict["aug11b"] = pkl.load(f)
with open("./sep1a/reactants_products_sep1a_filtered_cosmo_results.pkl", "rb") as f:
    cosmo_results_dict["sep1a"] = pkl.load(f)

def fill_column_cosmo(results_dict):
    global cosmo_result
    global df
    global props
    props = ['H (bar)', 'ln(gamma)', 'Pvap (bar)', 'Gsolv (kcal/mol)', 'Hsolv (kcal/mol)']
    for project, df in dfs_dict.items():
        cosmo_result = results_dict[project]

        solvent_names = cosmo_result["solvent_name"].unique()
        temps = cosmo_result["temp (K)"].unique()

        columns_dict = {f"{solvent_name}_{temp}_{prop}": [None for _ in df.index] for solvent_name in solvent_names for temp in temps for prop in props}

        out = Parallel(n_jobs=n_jobs, backend="multiprocessing", verbose=5)(delayed(create_column_cosmo)(row_ind) for row_ind in cosmo_result.index)

        for index, props_dict in out:
            for column_name, prop in props_dict.items():
                columns_dict[column_name][index] = prop

        for column_name, column in columns_dict.items():
            df[column_name] = column

def create_column_cosmo(row_ind):
    global cosmo_result
    global df
    global props
    props_dict = {}
    solvent_name = cosmo_result.loc[row_ind, "solvent_name"]
    temp = cosmo_result.loc[row_ind, "temp (K)"]
    mol_id = cosmo_result.loc[row_ind, "solute_name"]
    index = df.index[df["id"]==mol_id][0]
    for prop in props:
        props_dict[f"{solvent_name}_{temp}_{prop}"] = cosmo_result.loc[row_ind, prop]
    return index, props_dict

fill_column_cosmo(cosmo_results_dict)

df_merged = pd.concat(list(dfs_dict.values()), ignore_index=True)

with open("reactants_products_aug11b_sep1a_filtered_all_results_table.pkl", "wb") as f:
    pkl.dump(df_merged, f)