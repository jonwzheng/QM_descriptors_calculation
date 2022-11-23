import os
from rdkit import Chem
import subprocess
import shutil

from rdmc.mol import RDKitMol
from rdmc.forcefield import OpenBabelFF
from rdmc.ts import get_formed_and_broken_bonds

from lib.file_parser import xyz2com
from lib.semiempirical_calculation import run_xtb_opt
from lib.utils import mol2charge, mol2mult, mol2xyz

def reset_r_p_complex(rxn_smi, ts_xyz, ts_id, rdmc_path, g16_path, level_of_theory, n_procs, job_ram, subinputs_dir, suboutputs_dir, scratch_dir):
    current_dir = os.path.abspath(os.getcwd())

    r_complex_smi, p_complex_smi = rxn_smi.split(">>")
    r_complex = RDKitMol.FromSmiles(r_complex_smi)
    p_complex = RDKitMol.FromSmiles(p_complex_smi)
    ts_mol = RDKitMol.FromXYZ(ts_xyz, header=False)

    formed_bonds, broken_bonds = get_formed_and_broken_bonds(r_complex, p_complex)
    
    mol_id = f"{ts_id}_r"
    mol_scratch_dir = os.path.join(current_dir, mol_id)
    os.chdir(mol_scratch_dir)
    new_r_complex = reset_r_complex(ts_mol, r_complex, formed_bonds)
    xyz = mol2xyz(new_r_complex)
    charge = mol2charge(new_r_complex)
    mult = mol2mult(new_r_complex)
    run_xtb_opt(xyz, charge, mult, mol_id, rdmc_path, g16_path, level_of_theory, n_procs, job_ram, suboutputs_dir)
    os.chdir(current_dir)
    shutil.rmtree(mol_scratch_dir)

    mol_id = f"{ts_id}_p"
    mol_scratch_dir = os.path.join(current_dir, mol_id)
    new_p_complex = reset_p_complex(new_r_complex, p_complex, broken_bonds)
    xyz = mol2xyz(new_p_complex)
    charge = mol2charge(new_p_complex)
    mult = mol2mult(new_p_complex)
    run_xtb_opt(xyz, charge, mult, mol_id, rdmc_path, g16_path, level_of_theory, n_procs, job_ram, suboutputs_dir)
    os.chdir(current_dir)
    shutil.rmtree(mol_scratch_dir)

    os.remove(os.path.join(subinputs_dir, f"{ts_id}.tmp"))

def reset_r_complex(ts_mol, r_complex, formed_bonds):
    # copy current r_complex and set new positions
    new_r_complex = r_complex.Copy(quickCopy=True)
    new_r_complex.SetPositions(ts_mol.GetPositions())

    # setup first minimization with broken bond constraints
    obff = OpenBabelFF(force_field="uff")
    obff.setup(new_r_complex)
    ts_conf = ts_mol.GetConformer()
    current_distances = [ts_conf.GetBondLength(b) for b in formed_bonds]
    for b, d in zip(formed_bonds, current_distances):
        obff.add_distance_constraint(b, 1.5*d)
    obff.optimize(max_step=2000)

    # second minimization without constraints
    obff.constraints = None
    obff.optimize(max_step=2000)
    new_r_complex = obff.get_optimized_mol()

    # third optimization with MMFF94s
    obff = OpenBabelFF(force_field="mmff94s")
    obff.setup(new_r_complex)
    obff.optimize(max_step=2000)
    new_r_complex = obff.get_optimized_mol()
    
    return new_r_complex

def reset_p_complex(new_r_complex, p_complex, broken_bonds):
    # copy current p_complex and set new positions
    new_p_complex = p_complex.Copy(quickCopy=True)
    new_p_complex.SetPositions(new_r_complex.GetPositions())

    # setup first minimization with broken bond constraints
    obff = OpenBabelFF(force_field="uff")
    obff.setup(new_p_complex)
    r_conf = new_r_complex.GetConformer()
    current_distances = [r_conf.GetBondLength(b) for b in broken_bonds]
    for b, d in zip(broken_bonds, current_distances):
        obff.add_distance_constraint(b, 1.5*d)
    obff.optimize(max_step=2000)

    # second minimization without constraints
    obff.constraints = None
    obff.optimize(max_step=2000)
    new_p_complex = obff.get_optimized_mol()

    # third optimization with MMFF94s
    obff = OpenBabelFF(force_field="mmff94s")
    obff.setup(new_p_complex)
    obff.optimize(max_step=2000)
    new_p_complex = obff.get_optimized_mol()
    
    return new_p_complex
    