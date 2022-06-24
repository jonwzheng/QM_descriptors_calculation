from fileinput import filename
import os
import shutil
import subprocess

from rdkit import Chem
from .file_parser import mol2xyz


def dlpno_sp_calc(folder, sdf, orca_path, charge, mult, n_procs, job_ram):
    basename = os.path.basename(sdf)
    file_name = os.path.splitext(basename)[0]
    child_dir = os.path.abspath(os.path.join(folder, file_name))

    parent_dir = os.getcwd()

    os.chdir(child_dir)

    mol = Chem.SDMolSupplier(sdf, removeHs=False, sanitize=False)[0]
    xyz = mol2xyz(mol)
    coords = "\n".join(xyz.splitlines()[2:])
    
    script = generate_dlpno_sp_input(coords, charge, mult, job_ram, n_procs)

    infile = f"{file_name}.in"
    with open(infile, "w+") as f:
        f.write(script)

    #create working directory
    os.makedirs("scratch", exist_ok=True)
    os.chdir("scratch")
    shutil.copy(os.path.join(child_dir, infile), infile)

    #run jobs
    orca_command = os.path.join(orca_path, "orca")
    logfile = file_name + '.log'
    outfile = file_name + '.out'
    with open(outfile, 'w') as out:
        subprocess.run('{} {} > {}'.format(orca_command, infile, logfile), shell=True, stdout=out, stderr=out)

    os.chdir(child_dir)
    shutil.copy(os.path.join("scratch", logfile), logfile)
    shutil.copy(os.path.join("scratch", outfile), outfile)
    shutil.rmtree("scratch")

    # check for normal termination
    with open(logfile, "r") as f:
        lines = f.readlines()
        assert any(["ORCA TERMINATED NORMALLY" in line for line in lines])

    os.remove(sdf)
    os.chdir(parent_dir)

def generate_dlpno_sp_input(xyz_str: str,
                                charge: int,
                                multiplicity: int,
                                memory_mb: int,
                                cpu_threads: int,
                                ) -> str:
    """
    Modified from ACS
    """

    if multiplicity == 1:
        reference = 'rHF'
    elif multiplicity == 2:
        reference = 'uHF'
    else:
        raise NotImplementedError

    script = f"""!{reference} dlpno-ccsd(t) def2-tzvp def2-tzvp/c NormalPNO
!sp 

%maxcore {memory_mb}
%pal
nprocs {cpu_threads}
end

* xyz {charge} {multiplicity}
{xyz_str}
*"""
    return script

    