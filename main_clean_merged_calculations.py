import os
import shutil

submit_dir = os.getcwd()

# for suboutput_folder in os.listdir(os.path.join(submit_dir, "output", "FF_conf", "outputs")):
#     for job_folder in os.listdir(os.path.join(submit_dir, "output", "FF_conf", "outputs", suboutput_folder)):
#         mol_id = job_folder
#         if os.path.exists(os.path.join(submit_dir, "output", "FF_conf", "outputs", suboutput_folder, job_folder, mol_id+"_confs.sdf")):
#             shutil.move(os.path.join(submit_dir, "output", "FF_conf", "outputs", suboutput_folder, job_folder, mol_id+"_confs.sdf"), os.path.join(submit_dir, "output", "FF_conf", "outputs", suboutput_folder, mol_id+"_confs.sdf"))
#             print(os.path.join(submit_dir, "output", "FF_conf", "outputs", suboutput_folder, mol_id+"_confs.sdf"))
#             shutil.rmtree(os.path.join(submit_dir, "output", "FF_conf", "outputs", suboutput_folder, job_folder))

# for suboutput_folder in os.listdir(os.path.join(submit_dir, "output", "FF_conf", "outputs")):
#     if not any("_confs.sdf" in file for file in os.listdir(os.path.join(submit_dir, "output", "FF_conf", "outputs", suboutput_folder))):
#         print(os.path.join(submit_dir, "output", "FF_conf", "outputs", suboutput_folder))
        # shutil.rmtree(os.path.join(submit_dir, "output", "FF_conf", "outputs", suboutput_folder))

# for suboutput_folder in os.listdir(os.path.join(submit_dir, "output", "DFT_opt_freq", "outputs")):
#     for job_folder in os.listdir(os.path.join(submit_dir, "output", "DFT_opt_freq", "outputs", suboutput_folder)):
#         mol_id = job_folder
#         if os.path.exists(os.path.join(submit_dir, "output", "DFT_opt_freq", "outputs", suboutput_folder, job_folder, mol_id+".log")):
#             shutil.move(os.path.join(submit_dir, "output", "DFT_opt_freq", "outputs", suboutput_folder, job_folder, mol_id+".log"), os.path.join(submit_dir, "output", "DFT_opt_freq", "outputs", suboutput_folder, mol_id+".log"))
#             print(os.path.join(submit_dir, "output", "DFT_opt_freq", "outputs", suboutput_folder, mol_id+".log"))
#             shutil.rmtree(os.path.join(submit_dir, "output", "DFT_opt_freq", "outputs", suboutput_folder, job_folder))
