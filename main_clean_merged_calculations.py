import os
import shutil

submit_dir = os.getcwd()

for suboutput_folder in os.listdir(os.path.join(submit_dir, "output", "DFT_opt_freq", "outputs")):
    for job_folder in os.listdir(os.path.join(submit_dir, "output", "DFT_opt_freq", "outputs", suboutput_folder)):
        mol_id = job_folder
        if os.path.exists(os.path.join(submit_dir, "output", "DFT_opt_freq", "outputs", suboutput_folder, job_folder, mol_id+".log")):
            shutil.move(os.path.join(submit_dir, "output", "DFT_opt_freq", "outputs", suboutput_folder, job_folder, mol_id+".log"), os.path.join(submit_dir, "output", "DFT_opt_freq", "outputs", suboutput_folder, mol_id+".log"))
            print(os.path.join(submit_dir, "output", "DFT_opt_freq", "outputs", suboutput_folder, mol_id+".log"))
            shutil.rmtree(os.path.join(submit_dir, "output", "DFT_opt_freq", "outputs", suboutput_folder, job_folder))
