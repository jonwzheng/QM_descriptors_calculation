from argparse import ArgumentParser, Namespace
import os
import re
import shutil
import time
import yaml

import pandas as pd
import traceback

import rdkit.Chem as Chem

from lib import create_logger, done_jobs_record
from lib import csearch
from lib.xtb_optimization import xtb_optimization, xtb_status
from lib import dft_scf_qm_descriptor, dft_scf_opt
from lib import cosmo_calc
from lib import dlpno_sp_calc

parser = ArgumentParser()
parser.add_argument('--input_smiles', type=str, required=False,
                    help='input smiles included in a .csv file')
parser.add_argument('--output_folder', type=str, default='output',
                    help='output folder name')
parser.add_argument('--task_id', type=int, default=0,
                    help='task id for job arrays or LLsub')
parser.add_argument('--num_tasks', type=int, default=1,
                    help='Number of tasks for job arrays or LLsub')
# parser.add_argument('--output', type=str, default='QM_descriptors.pickle',
#                     help='output as a .pickle file')
# conformer searching
parser.add_argument('--conf_search_FF', type=str, default='GFNFF',
                    help='Force field that will be used for conformer search. Options are MMFF94s and GFNFF.')
parser.add_argument('--FF_conf_folder', type=str, default='FF_conf',
                    help='Folder name for FF searched conformers')
parser.add_argument('--nconf', type=int, default=800,
                    help='number of FF conformers')
parser.add_argument('-max_conf_try', type=int, default=2000,
                    help='maximum attempt for conformer generating, '
                         'this is useful for molecules with many chiral centers.')
parser.add_argument('-rmspre', type=float, required=False,
                        help='rms threshold pre optimization')
parser.add_argument('--rmspost', type=float, required=False, default=0.4,
                    help='rms threshold post FF minimization')
parser.add_argument('--E_cutoff_fraction', type=float, required=False, default=0.2,
                    help='energy window for FF minimization.')
parser.add_argument('--FF_threads', type=int, required=False, default=40,
                    help='number of process for the FF conformer searching')
parser.add_argument('--timeout', required=False, default=600,
                    help='time window for each FF conformer searching sub process')
parser.add_argument('--n_lowest_E_confs_to_save', type=int, default=10,
                    help='number of lowest energy conformers to save')

# xtb optimization and frequency calculation
parser.add_argument('--xtb_opt_freq_folder', type=str, default='XTB_opt_freq',
                    help='folder for XTB optimization')
parser.add_argument('--gaussian_xtb_opt_freq_theory', type=str, default='#opt=(calcall,maxcycle=128,noeig,nomicro)',
                    help='level of theory for the Gaussian-XTB calculation')
parser.add_argument('--gaussian_xtb_opt_n_procs', type=int, default=4,
                    help='number of process for Gaussian-XTB calculations')
parser.add_argument('--gaussian_xtb_opt_job_ram', type=int, default=3000,
                    help='amount of ram (MB) allocated for each Gaussian-XTB calculation')

# DFT optimization and frequency calculation
parser.add_argument('--DFT_opt_freq_folder', type=str, default='DFT_opt_freq',
                    help='folder for DFT optimization and frequency calculation',)
parser.add_argument('--DFT_opt_freq_theory', type=str, default='#p opt=(calcall,noeigentest,maxcycles=120) freq guess=mix wb97xd/def2svp scf=xqc iop(2/9=2000)',
                    help='level of theory for the DFT calculation')
parser.add_argument('--DFT_opt_freq_n_procs', type=int, default=4,
                    help='number of process for DFT calculations')
parser.add_argument('--DFT_opt_job_ram', type=int, default=3000,
                    help='amount of ram (MB) allocated for each DFT calculation')

# Turbomole and COSMO calculation
parser.add_argument('--COSMO_folder', type=str, default='COSMO_calc',
                    help='folder for COSMO calculation',)
parser.add_argument('--COSMO_temperature', type=float, required=True,
                    help='temperature used for COSMO calculation')
parser.add_argument('--COSMO_solvents', type=str, nargs="+", required=True,
                    help='solvents used for COSMO calculation')
parser.add_argument('--COSMO_solvent_solute_ratios', type=str, nargs="+", required=True,
                    help='solvent and solute ratios used for COSMO calculation')

# DLPNO single point calculation
parser.add_argument('--DLPNO_sp_folder', type=str, default='DLPNO_sp')
parser.add_argument('--DLPNO_sp_n_procs', type=int, default=4,
                    help='number of process for DLPNO calculations')
parser.add_argument('--DLPNO_sp_job_ram', type=int, default=3000,
                    help='amount of ram (MB) allocated for each DLPNO calculation')

# DFT QM descriptor calculation
parser.add_argument('--DFT_QM_descriptor_folder', type=str, default='DFT_QM_desc',
                    help='folder for DFT calculation')
parser.add_argument('--DFT_QM_descriptor_theory', type=str, default='b3lyp/def2svp',
                    help='level of theory for the DFT calculation')
parser.add_argument('--DFT_QM_descriptor_n_procs', type=int, default=4,
                    help='number of process for DFT calculations')
parser.add_argument('--DFT_QM_descriptor_job_ram', type=int, default=3000,
                    help='amount of ram (MB) allocated for each DFT calculation')

# Split number
parser.add_argument('--split', type=int, default=None,
                        help='split number for multi-part job')

# Job control
parser.add_argument('--only_DFT', action='store_true', help='only perform DFT related jobs')


# specify paths
parser.add_argument('--XTB_path', type=str, required=True,
                    help='path to installed XTB')
parser.add_argument('--G16_path', type=str, required=True,
                    help='path to installed Gaussian 16')
parser.add_argument('--RDMC_path', type=str, default=None,
                    help='path to RDMC to use xtb-gaussian script for xtb optimization calculation. If not provided, XTB will be used.')
parser.add_argument('--COSMOtherm_path', type=str, required=True,
                    help='path to COSMOthermo')
parser.add_argument('--COSMO_database_path', type=str, required=True,
                    help='path to COSMO_database')
parser.add_argument('--ORCA_path', type=str, required=True,
                    help='path to ORCA')

args = parser.parse_args()

XTB_PATH = args.XTB_path
G16_PATH = args.G16_path
RDMC_PATH = args.RDMC_path
COSMOTHERM_PATH = args.COSMOtherm_path
COSMO_DATABASE_PATH = args.COSMO_database_path
ORCA_PATH = args.ORCA_path

name = os.path.splitext(args.input_smiles)[0]
logger = create_logger(name=name, task_id=args.task_id)
project_dir = os.path.abspath(os.path.join(args.output_folder, f"{args.output_folder}_{args.task_id}"))

try:
    done_jobs_record.load(project_dir, args)
    logger.info("this is a restart job...")
    logger.info("loading completed job ids...")
except:
    logger.info("this is a new job...")
    logger.info("make output folder...")
    os.makedirs(args.output_folder, exist_ok=True)
    logger.info("making project folder...")
    os.makedirs(project_dir, exist_ok=True)

df = pd.read_csv(args.input_smiles, index_col=0)
df = df[args.task_id:len(df.index):args.num_tasks]
# create id to smile mapping
mol_id_to_smi_dict = dict(zip(df.id, df.smiles))
mol_id_to_charge_dict = dict()
mol_id_to_mult_dict = dict()
for k, v in mol_id_to_smi_dict.items():
    try:
        mol = Chem.MolFromSmiles(v)
    except Exception as e:
        logger.error(f'Cannot translate smi {v} to molecule for species {k}')

    try:
        charge = Chem.GetFormalCharge(mol)
        mol_id_to_charge_dict[k] = charge
    except Exception as e:
        logger.error(f'Cannot determine molecular charge for species {k} with smi {v}')

    num_radical_elec = 0
    for atom in mol.GetAtoms():
        num_radical_elec += atom.GetNumRadicalElectrons()
    mol_id_to_mult_dict[k] =  num_radical_elec + 1

# switch to project folder
logger.info("switching to project folder...")
os.chdir(project_dir)

# conformer searching

# if not args.only_DFT:
logger.info('starting FF conformer searching...')
supported_FFs = ["MMFF94s", "GFNFF"]
try:
    assert args.conf_search_FF in supported_FFs
except Exception as e:
    logger.error(f"{args.conf_search_FF} not in supported FFs.")
    raise
supp = (x for x in df[['id', 'smiles']].values if x[0] not in done_jobs_record.FF_conf)
done_jobs_record = csearch(supp, len(df), args, logger, done_jobs_record, project_dir)
conf_sdfs = [f"{mol_id}.sdf" for mol_id in done_jobs_record.FF_conf if mol_id not in done_jobs_record.XTB_opt_freq]
logger.info('='*80)
# xtb optimization

# if not args.only_DFT:
logger.info('starting GFN2-XTB structure optimization and frequency calculation for the lowest FF conformer...')
os.makedirs(args.xtb_opt_freq_folder, exist_ok=True)

if RDMC_PATH:
    logger.info("RDMC path provided. Using Gaussian-XTB to perform GFN2-XTB calculations.")
else:
    logger.info("RDMC path not provided. Using XTB to perform GFN2-XTB calculations.")

for conf_sdf in conf_sdfs:
    try:
        file_name = os.path.splitext(conf_sdf)[0]
        os.makedirs(os.path.join(args.xtb_opt_freq_folder, file_name), exist_ok=True)
        shutil.copyfile(os.path.join(args.FF_conf_folder, file_name, conf_sdf),
                        os.path.join(args.xtb_opt_freq_folder, file_name, conf_sdf))
        mol_id = file_name
        charge = mol_id_to_charge_dict[mol_id]
        mult = mol_id_to_mult_dict[mol_id]
        xtb_optimization(args.xtb_opt_freq_folder, conf_sdf, XTB_PATH, RDMC_PATH, G16_PATH, args.gaussian_xtb_opt_freq_theory, args.gaussian_xtb_opt_n_procs,
                                args.gaussian_xtb_opt_job_ram, charge, mult)
        logger.info(f'GFN2-XTB optimization and frequency calculation for {mol_id} completed')
        done_jobs_record.XTB_opt_freq.append(mol_id)
        done_jobs_record.save(project_dir, args)
    except Exception as e:
        logger.error('XTB optimization for {} failed'.format(os.path.splitext(conf_sdf)[0]))
        logger.error(traceback.format_exc())
        os.chdir(project_dir)
xtb_opt_sdfs = [f"{mol_id}_opt.sdf" for mol_id in done_jobs_record.XTB_opt_freq if mol_id not in done_jobs_record.DFT_opt_freq]
logger.info('GFN2-XTB optimization and frequency calculation finished.')

logger.info('='*80)
# else:
#     opt_sdfs = []
#     for mol_id, v in mol_id_to_smi_dict.items():
#         logger.info(f'checking xtb convergence for {mol_id}')
#         try:
#             opt_sdf = xtb_status(args.xtb_opt_freq_folder, mol_id)
#             opt_sdfs.append(opt_sdf)
#         except Exception as e:
#             logger.error('XTB optimization for {} failed: {}'.format(mol_id, e))

os.makedirs('yamls', exist_ok=True)

# G16 DFT calculation
# if not args.only_DFT:
#     os.makedirs(args.DFT_folder, exist_ok=True)
# else:
    # logger.info("Searching for optimized XTB files.")
    # opt_sdfs = []
    # for a_file in os.listdir(args.DFT_folder):
    #     if a_file.endswith(".sdf"):
    #         logger.info(f'Found file {a_file}')
    #         mol_id = a_file.split('_')[0]
    #         if mol_id in mol_id_to_smi_dict.keys():
    #             opt_sdfs.append(a_file)

logger.info('starting DFT optimization and frequency calculation for the XTB-optimized conformer...')
os.makedirs(args.DFT_opt_freq_folder, exist_ok=True)
for xtb_opt_sdf in xtb_opt_sdfs:
    try:
        file_name = os.path.splitext(xtb_opt_sdf)[0].split("_")[0]
        os.makedirs(os.path.join(args.DFT_opt_freq_folder, file_name), exist_ok=True)
        shutil.copyfile(os.path.join(args.xtb_opt_freq_folder, file_name, xtb_opt_sdf),
                        os.path.join(args.DFT_opt_freq_folder, file_name, file_name + ".sdf"))

        mol_id = file_name
        charge = mol_id_to_charge_dict[mol_id]
        mult = mol_id_to_mult_dict[mol_id]
        opt_sdf = dft_scf_opt(args.DFT_opt_freq_folder, file_name + ".sdf", G16_PATH, args.DFT_opt_freq_theory, args.DFT_opt_freq_n_procs,
                                logger, args.DFT_opt_job_ram, charge, mult)
        logger.info(f'DFT optimization and frequency calculation for {mol_id} completed')
        done_jobs_record.DFT_opt_freq.append(mol_id)
        done_jobs_record.save(project_dir, args)
    except Exception as e:
        logger.error('DFT optimization for {} failed'.format(os.path.splitext(conf_sdf)[0]))
        logger.error(traceback.format_exc())
        os.chdir(project_dir)
opt_sdfs = [f"{mol_id}_opt.sdf" for mol_id in done_jobs_record.DFT_opt_freq if mol_id not in done_jobs_record.COSMO]
logger.info('DFT optimization and frequency calculation finished.')
logger.info('='*80)

logger.info('starting Turbomole and COSMO calculation for the DFT-optimized conformer...')
os.makedirs(args.COSMO_folder, exist_ok=True)
for opt_sdf in opt_sdfs:
    try:
        file_name = os.path.splitext(opt_sdf)[0].split("_")[0]
        os.makedirs(os.path.join(args.COSMO_folder, file_name), exist_ok=True)
        shutil.copyfile(os.path.join(args.DFT_opt_freq_folder, file_name, opt_sdf),
                        os.path.join(args.COSMO_folder, file_name, file_name + ".sdf"))
        mol_id = file_name
        charge = mol_id_to_charge_dict[mol_id]
        mult = mol_id_to_mult_dict[mol_id]
        cosmo_calc(args.COSMO_folder, file_name + ".sdf", COSMOTHERM_PATH, COSMO_DATABASE_PATH, charge, mult, args.COSMO_solvents,
                   args.COSMO_solvent_solute_ratios, args.COSMO_temperature)
        done_jobs_record.COSMO.append(mol_id)
        done_jobs_record.save(project_dir, args)
        logger.info(f'COSMO calculation for {mol_id} completed')
    except:
        logger.error(f'Turbomole and COSMO calculation for {opt_sdf} failed.')
        logger.error(traceback.format_exc())
        os.chdir(project_dir)
logger.info('COSMO calculation finished.')
logger.info('='*80)

logger.info('starting DLPNO single point calculation for the DFT-optimized conformer...')
os.makedirs(args.DLPNO_sp_folder, exist_ok=True)
opt_sdfs = [f"{mol_id}_opt.sdf" for mol_id in done_jobs_record.COSMO if mol_id not in done_jobs_record.WFT_sp]
for opt_sdf in opt_sdfs:
    try:
        file_name = os.path.splitext(opt_sdf)[0].split("_")[0]
        os.makedirs(os.path.join(args.DLPNO_sp_folder, file_name), exist_ok=True)
        shutil.copyfile(os.path.join(args.DFT_opt_freq_folder, file_name, opt_sdf),
                        os.path.join(args.DLPNO_sp_folder, file_name, file_name + ".sdf"))
        mol_id = file_name
        charge = mol_id_to_charge_dict[mol_id]
        mult = mol_id_to_mult_dict[mol_id]
        dlpno_sp_calc(args.DLPNO_sp_folder, file_name + ".sdf", ORCA_PATH, charge, mult, args.DLPNO_sp_n_procs, args.DLPNO_sp_job_ram)
        done_jobs_record.WFT_sp.append(mol_id)
        done_jobs_record.save(project_dir, args)
        logger.info(f'DLPNO single point calculation for {mol_id} completed')
    except:
        logger.error(f'DLPNO single point calculation for {opt_sdf} failed.')
        logger.error(traceback.format_exc())
logger.info('DLPNO single point calculation finished.')


# # DFT QM descriptor calculation
# os.makedirs(args.DFT_QM_descriptor_folder, exist_ok=True)
# qm_descriptors = dict()
# for opt_sdf in opt_sdfs:
#     try:
#         mol_id = opt_sdf.split('_')[0]
#         charge = mol_id_to_charge_dict[mol_id]
#     except Exception as e:
#         logger.error(f'Cannot determine molecular charge for species {mol_id}')

#     # if not args.only_DFT:
#     try:
#         shutil.copyfile(os.path.join(args.xtb_opt_freq_folder, opt_sdf),
#                         os.path.join(args.DFT_QM_descriptor_folder, opt_sdf))
#         time.sleep(1)
#     except Exception as e:
#         logger.error(f'file IO error for {opt_sdf}')

# for opt_sdf in opt_sdfs:
#     _qm_descriptors = dict()
#     try:
#         mol_id = opt_sdf.split('_')[0]
#         charge = mol_id_to_charge_dict[mol_id]
#     except Exception as e:
#         logger.error(f'Cannot determine molecular charge for species {mol_id}')

#     try:
#         qm_descriptor = dft_scf_qm_descriptor(args.DFT_QM_descriptor_folder, opt_sdf, G16_PATH, args.DFT_QM_descriptor_theory, args.DFT_QM_descriptor_n_procs,
#                                 logger, args.DFT_QM_descriptor_job_ram, charge)
#     except Exception as e:
#         logger.error('Gaussian optimization for {} failed: {}'.format(os.path.splitext(opt_sdf)[0], e))

#     try:
#         mol_id = opt_sdf.split('_')[0]
#         smi = mol_id_to_smi_dict[mol_id]
#         qm_descriptors[mol_id] = (smi, qm_descriptor)
#         _qm_descriptors[mol_id] = (smi, qm_descriptor)
#         with open(f'yamls/{mol_id}_qm_descriptors.yaml', 'w') as output:
#             yaml.dump(_qm_descriptors, output)
#     except Exception as e:
#         logger.error(f'descriptor store error main.py line 143 - 144')

    

    

# if args.split is None:
#     with open('qm_descriptors.yaml', 'w') as output:
#         yaml.dump(qm_descriptors, output)
# else:
#     with open(f'qm_descriptors_{args.split}.yaml', 'w') as output:
#         yaml.dump(qm_descriptors, output)


