import os
import shutil

submit_dir = os.getcwd()
stages = ["FF_conf", "semiempirical_opt", "DFT_opt_freq", "COSMO_calc", "DLPNO_sp"]
splits = ["split_0", "split_1", "split_2", "test"]

os.makedirs(os.path.join(submit_dir, "output"), exist_ok=True)
for stage in stages:
    os.makedirs(os.path.join(submit_dir, "output", stage), exist_ok=True)
    os.makedirs(os.path.join(submit_dir, "output", stage, "outputs"), exist_ok=True)

for split in splits:
    if os.path.exists(os.path.join(submit_dir, split, "output")):
        for suboutput_folder in os.listdir(os.path.join(submit_dir, split, "output")):
            if "output_" in suboutput_folder:
                for stage in stages:
                    if os.path.exists(os.path.join(submit_dir, split, "output", suboutput_folder, stage)):
                        for job_folder in os.listdir(os.path.join(submit_dir, split, "output", suboutput_folder, stage)):
                            if "id" in job_folder:
                                ids = str(int(int(job_folder.split("id")[1])/1000))
                                os.makedirs(os.path.join(submit_dir, "output", stage, "outputs", f"outputs_{ids}"), exist_ok=True)
                                shutil.move(os.path.join(submit_dir, split, "output", suboutput_folder, stage, job_folder), os.path.join(submit_dir, "output", stage, "outputs", f"outputs_{ids}", job_folder))
                                print(os.path.join(submit_dir, "output", stage, "outputs", f"outputs_{ids}", job_folder))
