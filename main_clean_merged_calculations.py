import os
import shutil
import tarfile

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

# for suboutput_folder in os.listdir(os.path.join(submit_dir, "output", "semiempirical_opt", "outputs")):
#     for job_folder in os.listdir(os.path.join(submit_dir, "output", "semiempirical_opt", "outputs", suboutput_folder)):
#         mol_id = job_folder
#         if os.path.exists(os.path.join(submit_dir, "output", "semiempirical_opt", "outputs", suboutput_folder, job_folder, mol_id+".tar")):
#             shutil.move(os.path.join(submit_dir, "output", "semiempirical_opt", "outputs", suboutput_folder, job_folder, mol_id+".tar"), 
#                         os.path.join(submit_dir, "output", "semiempirical_opt", "outputs", suboutput_folder, mol_id+".tar"))
#             print(os.path.join(submit_dir, "output", "semiempirical_opt", "outputs", suboutput_folder, mol_id+".tar"))
#             shutil.rmtree(os.path.join(submit_dir, "output", "semiempirical_opt", "outputs", suboutput_folder, job_folder))

# for suboutput_folder in os.listdir(os.path.join(submit_dir, "output", "semiempirical_opt", "outputs")):
#     for job_folder in os.listdir(os.path.join(submit_dir, "output", "semiempirical_opt", "outputs", suboutput_folder)):
#         if os.path.isdir(os.path.join(submit_dir, "output", "semiempirical_opt", "outputs", suboutput_folder, job_folder)):
#             if any(".log" in log_file for log_file in os.listdir(os.path.join(submit_dir, "output", "semiempirical_opt", "outputs", suboutput_folder, job_folder))):
#                 mol_id = job_folder
#                 tar_file = os.path.join(submit_dir, "output", "semiempirical_opt", "outputs", suboutput_folder, f"{mol_id}.tar")
#                 tar = tarfile.open(tar_file, "w")
#                 for log_file in os.listdir(os.path.join(submit_dir, "output", "semiempirical_opt", "outputs", suboutput_folder, job_folder)):
#                     if ".log" in log_file:
#                         log_path = os.path.join(submit_dir, "output", "semiempirical_opt", "outputs", suboutput_folder, job_folder, log_file)
#                         tar.add(log_path)
#                 tar.close()
#                 print(tar_file)
#                 print(os.path.join(submit_dir, "output", "semiempirical_opt", "outputs", suboutput_folder, job_folder))
#                 shutil.rmtree(os.path.join(submit_dir, "output", "semiempirical_opt", "outputs", suboutput_folder, job_folder))

# for suboutput_folder in os.listdir(os.path.join(submit_dir, "output", "DFT_opt_freq", "outputs")):
#     for job_folder in os.listdir(os.path.join(submit_dir, "output", "DFT_opt_freq", "outputs", suboutput_folder)):
#         mol_id = job_folder
#         if os.path.exists(os.path.join(submit_dir, "output", "DFT_opt_freq", "outputs", suboutput_folder, job_folder, mol_id+".log")):
#             shutil.move(os.path.join(submit_dir, "output", "DFT_opt_freq", "outputs", suboutput_folder, job_folder, mol_id+".log"), os.path.join(submit_dir, "output", "DFT_opt_freq", "outputs", suboutput_folder, mol_id+".log"))
#             print(os.path.join(submit_dir, "output", "DFT_opt_freq", "outputs", suboutput_folder, mol_id+".log"))
#             shutil.rmtree(os.path.join(submit_dir, "output", "DFT_opt_freq", "outputs", suboutput_folder, job_folder))

# for suboutput_folder in os.listdir(os.path.join(submit_dir, "output", "DFT_opt_freq", "outputs")):
#     for job_folder in os.listdir(os.path.join(submit_dir, "output", "DFT_opt_freq", "outputs", suboutput_folder)):
#         if os.path.isdir(os.path.join(submit_dir, "output", "DFT_opt_freq", "outputs", suboutput_folder, job_folder)):
#             print(os.path.join(submit_dir, "output", "DFT_opt_freq", "outputs", suboutput_folder, job_folder))
#             shutil.rmtree(os.path.join(submit_dir, "output", "DFT_opt_freq", "outputs", suboutput_folder, job_folder))

# for suboutput_folder in os.listdir(os.path.join(submit_dir, "output", "COSMO_calc", "outputs")):
#     for job_folder in os.listdir(os.path.join(submit_dir, "output", "COSMO_calc", "outputs", suboutput_folder)):
#         mol_id = job_folder
#         if os.path.exists(os.path.join(submit_dir, "output", "COSMO_calc", "outputs", suboutput_folder, job_folder, mol_id+".tar")):
#             shutil.move(os.path.join(submit_dir, "output", "COSMO_calc", "outputs", suboutput_folder, job_folder, mol_id+".tar"), 
#                         os.path.join(submit_dir, "output", "COSMO_calc", "outputs", suboutput_folder, mol_id+".tar"))
#             print(os.path.join(submit_dir, "output", "COSMO_calc", "outputs", suboutput_folder, mol_id+".tar"))
#             shutil.rmtree(os.path.join(submit_dir, "output", "COSMO_calc", "outputs", suboutput_folder, job_folder))
            
for suboutput_folder in os.listdir(os.path.join(submit_dir, "output", "COSMO_calc", "outputs")):
    for job_folder in os.listdir(os.path.join(submit_dir, "output", "COSMO_calc", "outputs", suboutput_folder)):
        mol_id = job_folder
        if os.path.isdir(os.path.join(submit_dir, "output", "COSMO_calc", "outputs", suboutput_folder, job_folder)):
            print(os.path.join(submit_dir, "output", "COSMO_calc", "outputs", suboutput_folder, job_folder))
            # shutil.rmtree(os.path.join(submit_dir, "output", "COSMO_calc", "outputs", suboutput_folder, job_folder))