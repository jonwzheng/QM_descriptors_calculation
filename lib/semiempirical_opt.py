from rdkit import Chem
import os
import shutil
import subprocess
import traceback

import numpy as np

from .log_parser import XtbLog, G16Log
from .file_parser import mol2xyz, xyz2com, write_mol_to_sdf, write_mols_to_sdf


def semiempirical_opt(folder, sdf, xtb_path, rdmc_path, g16_path, level_of_theory, n_procs, job_ram, base_charge, mult, method, logger):
    basename = os.path.basename(sdf)
    file_name = os.path.splitext(basename)[0]

    parent_dir = os.getcwd()
    child_dir = os.path.abspath(os.path.join(folder, file_name))
    os.chdir(child_dir)

    mols = Chem.SDMolSupplier(sdf, removeHs=False, sanitize=False)
    os.remove(sdf)
    conf_ids_ens = []
    for conf_ind, mol in enumerate(mols):
        scratch_dir = f"{file_name}_{conf_ind}"
        os.makedirs(scratch_dir, exist_ok=True)
        os.chdir(scratch_dir)

        xyz = mol2xyz(mol)

        g16_command = os.path.join(g16_path, 'g16')

        if method == "GFN2-XTB":
            head = '%nprocshared={}\n%mem={}mb\n{}\nexternal=\"{}/rdmc/external/xtb_tools/xtb_gaussian.pl --gfn 2 -P\"\n'.format(n_procs, job_ram, level_of_theory, rdmc_path)
        else:
            head = '%nprocshared={}\n%mem={}mb\n{} {}\n'.format(n_procs, job_ram, level_of_theory, method)

        comfile = f"{file_name}_{conf_ind}.gjf"
        xyz2com(xyz, head=head, comfile=comfile, charge=base_charge, mult=mult, footer='\n')

        logfile = f"{file_name}_{conf_ind}.log"
        outfile = f"{file_name}_{conf_ind}.out"

        with open(outfile, 'w') as out:
            subprocess.run('{} < {} >> {}'.format(g16_command, comfile, logfile), shell=True, stdout=out, stderr=out)

        log = G16Log(logfile)
        if log.termination and np.min(log.har_frequencies) > 0:
            conf_ids_ens.append((conf_ind, log.E))
            conf = mol.GetConformer()
            for i in range(mol.GetNumAtoms()):
                conf.SetAtomPosition(i, log.Coords[i,:])
        else:
            logger.error(f'optimization of conformer {conf_ind} for {file_name} failed.')
        os.chdir(child_dir)
        shutil.rmtree(scratch_dir)

    write_mols_to_sdf(mols, f'{file_name}_confs.sdf')
    conf_ids_ens.sort(key=lambda x: x[1])
    write_mol_to_sdf(mols[conf_ids_ens[0][0]], f'{file_name}_opt.sdf')
    
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