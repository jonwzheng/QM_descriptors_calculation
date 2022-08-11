from fileinput import filename
import os
import shutil
import subprocess

from rdkit import Chem
from .file_parser import mol2xyz


def dlpno_sp_calc(mol_id, orca_path, charge, mult, n_procs, job_ram, xyz_DFT_opt):
    sdf = mol_id + '.sdf'
    mol_dir = os.getcwd()

    if xyz_DFT_opt:
        coords = xyz_DFT_opt[mol_id]
    else:
        mol = Chem.SDMolSupplier(sdf, removeHs=False, sanitize=False)[0]
        xyz = mol2xyz(mol)
        coords = "\n".join(xyz.splitlines()[2:])
    
    script = generate_dlpno_sp_input(coords, charge, mult, job_ram, n_procs)

    infile = f"{mol_id}.in"
    with open(infile, "w+") as f:
        f.write(script)

    #create working directory
    os.makedirs("scratch", exist_ok=True)
    os.chdir("scratch")
    shutil.copy(os.path.join(mol_dir, infile), infile)

    #run jobs
    orca_command = os.path.join(orca_path, "orca")
    logfile = mol_id + '.log'
    outfile = mol_id + '.out'
    with open(outfile, 'w') as out:
        subprocess.run('{} {} > {}'.format(orca_command, infile, logfile), shell=True, stdout=out, stderr=out)

    os.chdir(mol_dir)
    shutil.copy(os.path.join("scratch", logfile), logfile)
    shutil.copy(os.path.join("scratch", outfile), outfile)
    shutil.rmtree("scratch")

    # check for normal termination
    with open(logfile, "r") as f:
        lines = f.readlines()
        assert any(["ORCA TERMINATED NORMALLY" in line for line in lines])

    os.remove(sdf)

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

    