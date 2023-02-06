import os
import sys
import pandas as pd

# make arkane species and thermo input, and run arkane

task_id = int(sys.argv[1])
num_tasks = int(sys.argv[2])

input_df_paths = {
    "aug11b": "../calculations/aug11b/inputs/reactants_products_aug11b_inputs.csv",
    "sep1a_filtered": "../calculations/sep1a_filtered/inputs/reactants_products_sep1a_filtered_inputs.csv"
}

# make output directories
output_dir = "output"
os.makedirs(output_dir, exist_ok=True)

for project_name, input_df_path in input_df_paths.items():

    # make project directories
    project_dir = os.path.join(output_dir, project_name)
    os.makedirs(project_dir, exist_ok=True)

    inputs_dir = os.path.join(project_dir, "inputs")
    os.makedirs(inputs_dir, exist_ok=True)
    outputs_dir = os.path.join(project_dir, "outputs")
    os.makedirs(outputs_dir, exist_ok=True)


    input_df = pd.read_csv(input_df_path)
    mol_ids = list(input_df.id)
    mol_smis = list(input_df.smiles)
    mol_ids_mol_smis = list(zip(mol_ids, mol_smis))

    # make helper input files
    for mol_id, mol_smi in mol_ids_mol_smis[task_id::num_tasks]:

        ids = str(int(int(mol_id.split("id")[1])/1000))

        subinputs_dir = os.path.join(inputs_dir, f"inputs_{ids}")
        os.makedirs(subinputs_dir, exist_ok=True)
        suboutputs_dir = os.path.join(outputs_dir, f"outputs_{ids}")
        os.makedirs(suboutputs_dir, exist_ok=True)

        suboutputs_mol_dir = os.path.join(suboutputs_dir, mol_id)
        os.makedirs(suboutputs_mol_dir, exist_ok=True)
        arkane_output_file = os.path.join(suboutputs_mol_dir, "output.py")
        
        if not os.path.exists(arkane_output_file):
            
            os.makedirs(subinputs_dir, exist_ok=True)
            input_file = os.path.join(subinputs_dir, f"{mol_id}.in")

            if not os.path.exists(input_file):

                with open(input_file, "w+") as f:
                    f.write("")

    # run arkane
    for subinputs_folder in os.listdir(inputs_dir):
        subinputs_dir = os.path.join(inputs_dir, subinputs_folder)
        
        for input_file in os.listdir(subinputs_dir):
            if input_file.endswith(".in"):
                
                mol_id = input_file.split(".")[0]
                try:
                    os.rename()


