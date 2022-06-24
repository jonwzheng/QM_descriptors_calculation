from rdkit import Chem
import os
import shutil
import subprocess
import traceback

import numpy as np

from .log_parser import XtbLog, G16Log
from .file_parser import mol2xyz, xyz2com, write_mol_to_sdf


def xtb_optimization(folder, sdf, xtb_path, rdmc_path, g16_path, level_of_theory, n_procs, job_ram, base_charge, mult):
    basename = os.path.basename(sdf)
    file_name = os.path.splitext(basename)[0]
    scratch_dir = os.path.join(folder, file_name)

    parent_dir = os.getcwd()

    os.chdir(scratch_dir)

    if rdmc_path:
        mol = Chem.SDMolSupplier(sdf, removeHs=False, sanitize=False)[0]
        xyz = mol2xyz(mol)

        g16_command = os.path.join(g16_path, 'g16')
        head = '%nprocshared={}\n%mem={}mb\n{}\nexternal=\"{}/rdmc/external/xtb_tools/xtb_gaussian.pl --gfn 2 -P\"\n'.format(n_procs, job_ram, level_of_theory, rdmc_path)

        comfile = file_name + '.gjf'
        xyz2com(xyz, head=head, comfile=comfile, charge=base_charge, mult=mult, footer='\n')

        logfile = file_name + '.log'
        outfile = file_name + '.out'

        with open(outfile, 'w') as out:
            subprocess.run('{} < {} >> {}'.format(g16_command, comfile, logfile), shell=True, stdout=out, stderr=out)
            os.remove(sdf)

        log = G16Log(logfile)
        conf = mol.GetConformer()
        for i in range(mol.GetNumAtoms()):
            conf.SetAtomPosition(i, log.Coords[i,:])
        write_mol_to_sdf(mol, f'{file_name}_opt.sdf')
        
    else:
        xtb_command = os.path.join(xtb_path, 'xtb')
        with open(os.path.join('{}_opt.log'.format(file_name)), 'w') as out:
            subprocess.call([xtb_command, os.path.join('{}.sdf'.format(file_name)), '-opt'],
                            stdout=out, stderr=out)
            shutil.move('xtbopt.sdf', os.path.join('{}_opt.sdf'.format(file_name)))
            os.remove(sdf)

        with open(os.path.join(file_name + '_freq.log'), 'w') as out:
            subprocess.call([xtb_command, os.path.join('{}_opt.sdf'.format(file_name)), '-ohess'], stdout=out,
                            stderr=out)

        log = XtbLog(os.path.join('{}_freq.log'.format(file_name)))
    if log.termination:
        if rdmc_path:
            peaks = log.har_frequencies
        else:
            peaks = log.wavenum
        if np.min(peaks) < 0:
            raise RuntimeError('imaginary frequency found for {}'.format(file_name))
    else:
        raise RuntimeError('xtb optimization did not finish for {}'.format(file_name))
    
    os.chdir(parent_dir)

def xtb_status(folder, molid):

    try:
        log = XtbLog(os.path.join(folder,'{}_freq.log'.format(molid)))
    except:
        raise RuntimeError(f'xtb log file not found for {molid}')

    if log.termination:
        peaks = log.wavenum
        if np.min(peaks) < 0:
            raise RuntimeError('imaginary frequency found for {}'.format(molid))
        else:
            return '{}_opt.sdf'.format(molid)
    else:
        raise RuntimeError('xtb optimization did not finish for {}'.format(molid))