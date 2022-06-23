from fileinput import filename
import os
import shutil
import subprocess

from rdkit import Chem
from .file_parser import mol2xyz


def cosmo_calc(folder, sdf, cosmotherm_path, cosmo_database_path, charge, mult, solvents, ratios, temperature):
    basename = os.path.basename(sdf)
    file_name = os.path.splitext(basename)[0]
    parent_dir = os.getcwd()
    
    child_dir = os.path.abspath(os.path.join(folder, file_name))
    os.chdir(child_dir)

    mol = Chem.SDMolSupplier(sdf, removeHs=False, sanitize=False)[0]
    xyz = mol2xyz(mol)

    # prepare for turbomole calculation
    os.makedirs("xyz", exist_ok=True)
    xyz_file_name = f'{file_name}.xyz'
    with open(os.path.join("xyz", xyz_file_name), "w+") as f:
        f.write(xyz)

    txtfile = f'{file_name}.txt'
    with open(txtfile, "w+") as f:
        f.write(f"{file_name} {charge} {mult}")

    #create and move to working directory
    os.makedirs("scratch", exist_ok=True)
    os.chdir("scratch")
    shutil.copy(os.path.join(child_dir, txtfile), txtfile)
    shutil.copytree(os.path.join(child_dir, "xyz"), "xyz")

    #run the job
    logfile = file_name + '.log'
    outfile = file_name + '.out'
    with open(outfile, 'w') as out:
        subprocess.run(f'calculate -l {txtfile} -m BP-TZVPD-FINE-COSMO-SP -f xyz -din xyz > {logfile}', shell=True, stdout=out, stderr=out)
        subprocess.run(f'calculate -l {txtfile} -m BP-TZVPD-GAS-SP -f xyz -din xyz > {logfile}', shell=True, stdout=out, stderr=out)

    #move files back
    shutil.copy(logfile, os.path.join(child_dir, logfile))
    shutil.copy(outfile, os.path.join(child_dir, outfile))

    for file in os.listdir("CosmofilesBP-TZVPD-FINE-COSMO-SP"):
        if file.endswith("cosmo"):
            shutil.copy(os.path.join("CosmofilesBP-TZVPD-FINE-COSMO-SP",file), file)
            shutil.copy(os.path.join("CosmofilesBP-TZVPD-FINE-COSMO-SP",file), os.path.join(child_dir, file))
    for file in os.listdir("EnergyfilesBP-TZVPD-FINE-COSMO-SP"):
        if file.endswith("energy"):
            shutil.copy(os.path.join("EnergyfilesBP-TZVPD-FINE-COSMO-SP", file), file)
            shutil.copy(os.path.join("EnergyfilesBP-TZVPD-FINE-COSMO-SP", file), os.path.join(child_dir, file))
    
    # prepare for cosmo calculation
    script = generate_cosmo_input(file_name, solvents, ratios, temperature, cosmotherm_path, cosmo_database_path)

    inpfile = f'{file_name}.inp'
    with open(inpfile, "w+") as f:
        f.write(script)

    cosmo_command = os.path.join(cosmotherm_path, "COSMOtherm", "BIN-LINUX", "cosmotherm")
    outfile = file_name + '.out'
    tabfile = file_name + '.tab'
    with open(outfile, 'w') as out:
        subprocess.run(f'{cosmo_command} {inpfile}', shell=True, stdout=out, stderr=out)

    #move files back
    shutil.copy(inpfile, os.path.join(child_dir, inpfile))
    shutil.copy(outfile, os.path.join(child_dir, outfile))
    shutil.copy(tabfile, os.path.join(child_dir, tabfile))

    os.chdir(child_dir)

    #remove working directory
    shutil.rmtree("scratch")
    os.remove(sdf)
    os.chdir(parent_dir)
    
def generate_cosmo_input(name, solvents, ratios, temperature, cosmotherm_path, cosmo_database_path):
    """
    Modified from ACS
    """

    script = f"""ctd = BP_TZVPD_FINE_20.ctd cdir = "{cosmotherm_path}/COSMOtherm/CTDATA-FILES" ldir = "{cosmotherm_path}/licensefiles"
unit notempty wtln ehfile
!! Title !!
"""

    for solvent in solvents:
        script += f"""f = "{solvent}_c0.cosmo" fdir="{cosmo_database_path}/COSMObase2020/BP-TZVPD-FINE/{solvent[0]}" VPfile
"""

    script += f"""f = "{name}.cosmo" fdir="." VPfile
henry  xh={{ {" ".join(ratios)} }}  tc={temperature} GSOLV # Automatic Henry Law coefficient Calculation
"""
    return script