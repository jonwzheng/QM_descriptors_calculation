from argparse import ArgumentParser
import os
import shutil

import pandas as pd
import traceback

import rdkit.Chem as Chem

from lib import create_logger, DoneJobsRecord
from lib import csearch
from lib import semiempirical_opt
from lib import dft_scf_qm_descriptor, dft_scf_opt, dft_scf_sp, save_dft_sp_results
from lib import cosmo_calc, save_cosmo_results
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
parser.add_argument('--max_n_conf', type=int, default=800,
                    help='maximum number of FF conformers. nc = 3**n_rotatable_bonds, n_conf = nc if nc < max_n_conf else max_n_conf')
parser.add_argument('-max_conf_try', type=int, default=2000,
                    help='maximum attempt for conformer generating, '
                         'this is useful for molecules with many chiral centers.')
parser.add_argument('-rmspre', type=float, required=False, default=0.1,
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

# semiempirical optimization and frequency calculation
parser.add_argument('--semiempirical_opt_folder', type=str, default='semiempirical_opt',
                    help='folder for semiempirical optimization')
parser.add_argument('--semiempirical_method', type=str, default='GFN2-XTB',
                    help='method used for semiempirical optimization. Options are GFN2-XTB, am1, and pm7.')
parser.add_argument('--gaussian_semiempirical_opt_theory', type=str, default='#opt=(calcall,maxcycle=128,noeig,nomicro)',
                    help='level of theory for the Gaussian semiempirical calculation')
parser.add_argument('--gaussian_semiempirical_opt_n_procs', type=int, default=4,
                    help='number of process for Gaussian semiempirical calculations')
parser.add_argument('--gaussian_semiempirical_opt_job_ram', type=int, default=3000,
                    help='amount of ram (MB) allocated for each Gaussian semiempirical calculation')

# DFT optimization and frequency calculation
parser.add_argument('--DFT_opt_freq_folder', type=str, default='DFT_opt_freq',
                    help='folder for DFT optimization and frequency calculation',)
parser.add_argument('--DFT_opt_freq_theory', type=str, default='#p opt=(calcall,noeigentest,maxcycles=120) freq guess=mix wb97xd/def2svp scf=xqc iop(2/9=2000)',
                    help='level of theory for the DFT calculation')
parser.add_argument('--DFT_opt_freq_n_procs', type=int, default=4,
                    help='number of process for DFT calculations')
parser.add_argument('--DFT_opt_freq_job_ram', type=int, default=3000,
                    help='amount of ram (MB) allocated for each DFT calculation')

# Turbomole and COSMO calculation
parser.add_argument('--COSMO_folder', type=str, default='COSMO_calc',
                    help='folder for COSMO calculation',)
parser.add_argument('--COSMO_temperatures', type=str, nargs="+", required=False, default=['297.15', '298.15', '299.15'],
                    help='temperatures used for COSMO calculation')
parser.add_argument('--COSMO_input_pure_solvents', type=str, required=False, default='common_solvent_list_final.csv',
                    help='input file containing pure solvents used for COSMO calculation.')

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

# test
parser.add_argument('--is_test', type=bool, default=False,
                    help='whether this is to test different semiempirical methods')
# DFT single point calculation for test
parser.add_argument('--DFT_sp_folder', type=str, default='DFT_sp',
                    help='folder for DFT optimization and frequency calculation',)
parser.add_argument('--DFT_sp_theory', type=str, default='#p guess=mix wb97xd/def2svp scf=xqc iop(2/9=2000)',
                    help='level of theory for the DFT calculation')
parser.add_argument('--DFT_sp_n_procs', type=int, default=4,
                    help='number of process for DFT calculations')
parser.add_argument('--DFT_sp_job_ram', type=int, default=3000,
                    help='amount of ram (MB) allocated for each DFT calculation')

# specify paths
parser.add_argument('--XTB_path', type=str, required=True,
                    help='path to installed XTB')
parser.add_argument('--G16_path', type=str, required=True,
                    help='path to installed Gaussian 16')
parser.add_argument('--RDMC_path', type=str, required=True,
                    help='path to RDMC to use xtb-gaussian script for xtb optimization calculation.')
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
submit_dir = os.path.abspath(os.getcwd())
project_dir = os.path.abspath(os.path.join(args.output_folder, f"{args.output_folder}_{args.task_id}"))

done_jobs_record = DoneJobsRecord()

try:
    done_jobs_record.load(project_dir, args.task_id)
    logger.info("this is a restart job...")
    logger.info("loading completed job ids...")
except:
    logger.info("this is a new job...")
    logger.info("make output folder...")
    os.makedirs(args.output_folder, exist_ok=True)
    logger.info("making project folder...")
    os.makedirs(project_dir, exist_ok=True)

df = pd.read_csv(args.input_smiles, index_col=0)
assert len(df['id']) == len(set(df['id'])), "ids must be unique"
df.sort_values(by='smiles', key=lambda x: x.str.len(), inplace=True) #sort by length of smiles to help even out the workload of each task
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
logger.info('starting FF conformer searching...')
supported_FFs = ["MMFF94s", "GFNFF"]
try:
    assert args.conf_search_FF in supported_FFs
except Exception as e:
    logger.error(f"{args.conf_search_FF} not in supported FFs.")
    raise
supp = (x for x in df[['id', 'smiles']].values if x[0] not in done_jobs_record.FF_conf)
conf_ids = [x[0] for x in df[['id', 'smiles']].values if x[0] not in done_jobs_record.FF_conf]
if conf_ids:
    conf_ids_str = ','.join(conf_ids)
    logger.info(f'FF conformer searching for: {conf_ids_str}')
    done_jobs_record = csearch(supp, len(conf_ids), args, logger, done_jobs_record, project_dir)

logger.info('FF conformer searching completed.')
logger.info('='*80)

if args.is_test:
    semiempirical_methods = ["GFN2-XTB", "am1", "pm7"]
    conf_sdfs = [f"{mol_id}_confs.sdf" for mol_id in done_jobs_record.FF_conf if len(done_jobs_record.test_semiempirical_opt.get(mol_id, []))!= len(semiempirical_methods)]
    logger.info(f'starting geometry optimization for lowest energy FF-optimized conformers using different semiempirical methods...')
    os.makedirs(args.semiempirical_opt_folder, exist_ok=True)

    for conf_sdf in conf_sdfs:
        mol_id = os.path.splitext(conf_sdf)[0].split("_")[0]
        os.makedirs(os.path.join(args.semiempirical_opt_folder, mol_id), exist_ok=True)
        charge = mol_id_to_charge_dict[mol_id]
        mult = mol_id_to_mult_dict[mol_id]

        for semiempirical_method in semiempirical_methods:
            if semiempirical_method not in done_jobs_record.test_semiempirical_opt.get(mol_id, []):
                os.makedirs(os.path.join(args.semiempirical_opt_folder, mol_id, semiempirical_method), exist_ok=True)
                shutil.copyfile(os.path.join(args.FF_conf_folder, mol_id, conf_sdf),
                                os.path.join(args.semiempirical_opt_folder, mol_id, semiempirical_method, mol_id + ".sdf"))

                os.chdir(os.path.join(args.semiempirical_opt_folder, mol_id, semiempirical_method))
                logger.info(f'starting {semiempirical_method} semiempirical geometry optimization calculation for {mol_id}...') 
                try:
                    semiempirical_opt(mol_id, XTB_PATH, RDMC_PATH, G16_PATH, args.gaussian_semiempirical_opt_theory, args.gaussian_semiempirical_opt_n_procs,
                                                args.gaussian_semiempirical_opt_job_ram, charge, mult, semiempirical_method, logger)
                    done_jobs = done_jobs_record.test_semiempirical_opt.get(mol_id, [])
                    done_jobs.append(semiempirical_method)
                    done_jobs_record.test_semiempirical_opt[mol_id] = done_jobs
                    done_jobs_record.save(project_dir, args.task_id)
                    logger.info(f'all {semiempirical_method} semiempirical geometry optimization calculations for {mol_id} completed')
                except:
                    logger.error(f'{semiempirical_method} semiempirical geometry optimization calculation for {mol_id} failed')
                    logger.error(traceback.format_exc())
                os.chdir(project_dir)

    semi_opt_sdfs = [f"{mol_id}_opt.sdf" for mol_id in done_jobs_record.test_semiempirical_opt if len(done_jobs_record.test_DFT_sp.get(mol_id, []))!= len(semiempirical_methods)]
    logger.info(f'semiempirical geometry optimization calculations finished.')
    logger.info('='*80)

    logger.info('starting DFT single point calculations for the lowest energy semiempirical-optimized conformer...')
    os.makedirs(os.path.join(args.DFT_sp_folder), exist_ok=True)
    for semi_opt_sdf in semi_opt_sdfs:
        mol_id = os.path.splitext(semi_opt_sdf)[0].split("_")[0]
        os.makedirs(os.path.join(args.DFT_sp_folder, mol_id), exist_ok=True)
        charge = mol_id_to_charge_dict[mol_id]
        mult = mol_id_to_mult_dict[mol_id]

        for semiempirical_method in semiempirical_methods:
            if semiempirical_method in done_jobs_record.test_semiempirical_opt.get(mol_id, []) and semiempirical_method not in done_jobs_record.test_DFT_sp.get(mol_id, []):
                os.makedirs(os.path.join(args.DFT_sp_folder, mol_id, semiempirical_method), exist_ok=True)
                shutil.copyfile(os.path.join(args.semiempirical_opt_folder, mol_id, semiempirical_method, semi_opt_sdf),
                                os.path.join(args.DFT_sp_folder, mol_id, semiempirical_method, mol_id + ".sdf"))

                os.chdir(os.path.join(args.DFT_sp_folder, mol_id, semiempirical_method))
                try:
                    dft_scf_sp(mol_id, G16_PATH, args.DFT_sp_theory, args.DFT_sp_n_procs, logger, args.DFT_sp_job_ram, charge, mult)
                    done_jobs = done_jobs_record.test_DFT_sp.get(mol_id, [])
                    done_jobs.append(semiempirical_method)
                    done_jobs_record.test_DFT_sp[mol_id] = done_jobs
                    done_jobs_record.save(project_dir, args.task_id)
                    logger.info(f'DFT single point calculation for {semiempirical_method} optimized {mol_id} completed')
                except:
                    logger.error(f'DFT single point calculation for {semiempirical_method} optimized {mol_id} failed')
                    logger.error(traceback.format_exc())
                os.chdir(project_dir)
    logger.info('DFT single point calculations finished.')
    logger.info('='*80)
    logger.info('All calculations completed.')
    logger.info('Extracting DFT single point calculation results for comparison...')
    save_dft_sp_results(args.DFT_sp_folder, done_jobs_record, args.task_id, mol_id_to_smi_dict, semiempirical_methods)
    logger.info('Extrcating DFT single point calculation results completed.')

else:
    conf_sdfs = [f"{mol_id}_confs.sdf" for mol_id in done_jobs_record.FF_conf if mol_id not in done_jobs_record.semiempirical_opt]

    # semiempirical optimization
    logger.info(f'starting semiempirical geometry optimization for lowest energy FF-optimized conformers...')
    os.makedirs(args.semiempirical_opt_folder, exist_ok=True)

    for conf_sdf in conf_sdfs:
        mol_id = os.path.splitext(conf_sdf)[0].split("_")[0]
        os.makedirs(os.path.join(args.semiempirical_opt_folder, mol_id), exist_ok=True)
        shutil.copyfile(os.path.join(args.FF_conf_folder, mol_id, conf_sdf),
                        os.path.join(args.semiempirical_opt_folder, mol_id, mol_id + ".sdf"))
        charge = mol_id_to_charge_dict[mol_id]
        mult = mol_id_to_mult_dict[mol_id]
        os.chdir(os.path.join(args.semiempirical_opt_folder, mol_id))
        try:
            semiempirical_opt(mol_id, XTB_PATH, RDMC_PATH, G16_PATH, args.gaussian_semiempirical_opt_theory, args.gaussian_semiempirical_opt_n_procs,
                              args.gaussian_semiempirical_opt_job_ram, charge, mult, args.semiempirical_method, logger)
            done_jobs_record.semiempirical_opt.append(mol_id)
            done_jobs_record.save(project_dir, args.task_id)
            logger.info(f'semiempirical optimization for {mol_id} completed')
        except Exception as e:
            logger.error(f'semiempirical optimization for {mol_id} failed')
            logger.error(traceback.format_exc())
        os.chdir(project_dir)
    semi_opt_sdfs = [f"{mol_id}_opt.sdf" for mol_id in done_jobs_record.semiempirical_opt if mol_id not in done_jobs_record.DFT_opt_freq]
    logger.info('semiempirical optimization finished.')
    logger.info('='*80)

    logger.info('starting DFT optimization and frequency calculation for the lowest energy semiempirical-optimized conformer...')
    os.makedirs(args.DFT_opt_freq_folder, exist_ok=True)
    for semi_opt_sdf in semi_opt_sdfs:
        mol_id = os.path.splitext(semi_opt_sdf)[0].split("_")[0]
        os.makedirs(os.path.join(args.DFT_opt_freq_folder, mol_id), exist_ok=True)
        shutil.copyfile(os.path.join(args.semiempirical_opt_folder, mol_id, semi_opt_sdf),
                        os.path.join(args.DFT_opt_freq_folder, mol_id, mol_id + ".sdf"))

        charge = mol_id_to_charge_dict[mol_id]
        mult = mol_id_to_mult_dict[mol_id]
        os.chdir(os.path.join(args.DFT_opt_freq_folder, mol_id))
        try:
            dft_scf_opt(mol_id, G16_PATH, args.DFT_opt_freq_theory, args.DFT_opt_freq_n_procs,
                        logger, args.DFT_opt_freq_job_ram, charge, mult)
            done_jobs_record.DFT_opt_freq.append(mol_id)
            done_jobs_record.save(project_dir, args.task_id)
            logger.info(f'DFT optimization and frequency calculation for {mol_id} completed')
        except Exception as e:
            logger.error(f'DFT optimization and frequency calculation for {mol_id} failed')
            logger.error(traceback.format_exc())
        os.chdir(project_dir)
    logger.info('DFT optimization and frequency calculation finished.')
    logger.info('='*80)

    logger.info('starting Turbomole and COSMO calculation for the DFT-optimized conformer...')
    os.makedirs(args.COSMO_folder, exist_ok=True)
    logger.info('load solvent file...')
    df_pure = pd.read_csv(os.path.join(submit_dir,args.COSMO_input_pure_solvents))
    df_pure = df_pure.reset_index()
    opt_sdfs = [f"{mol_id}_opt.sdf" for mol_id in done_jobs_record.DFT_opt_freq if len(done_jobs_record.COSMO.get(mol_id, [])) < len(df_pure.index)]

    for opt_sdf in opt_sdfs:
        mol_id = os.path.splitext(opt_sdf)[0].split("_")[0]
        os.makedirs(os.path.join(args.COSMO_folder, mol_id), exist_ok=True)
        shutil.copyfile(os.path.join(args.DFT_opt_freq_folder, mol_id, opt_sdf),
                        os.path.join(args.COSMO_folder, mol_id, mol_id + ".sdf"))
        charge = mol_id_to_charge_dict[mol_id]
        mult = mol_id_to_mult_dict[mol_id]
        os.chdir(os.path.join(args.COSMO_folder, mol_id))
        try:
            cosmo_calc(mol_id, COSMOTHERM_PATH, COSMO_DATABASE_PATH, charge, mult, args.COSMO_temperatures, df_pure, done_jobs_record, project_dir, args.task_id)
            done_jobs = done_jobs_record.COSMO.get(mol_id, [])
            done_jobs.append(mol_id)
            done_jobs_record.COSMO[mol_id] = done_jobs
            logger.info(f'COSMO calculation for {mol_id} completed')
        except:
            logger.error(f'Turbomole and COSMO calculation for {opt_sdf} failed.')
            logger.error(traceback.format_exc())
        os.chdir(project_dir)

    logger.info('COSMO calculation finished.')
    logger.info('Extracting COSMO results...')
    try:
        save_cosmo_results(args.COSMO_folder, done_jobs_record, args.task_id)
        logger.error('Extractomg COSMO results completed.')
    except:
        logger.error('Extractomg COSMO results failed.')
        logger.error(traceback.format_exc())

    logger.info('='*80)

    logger.info('starting DLPNO single point calculation for the DFT-optimized conformer...')
    os.makedirs(args.DLPNO_sp_folder, exist_ok=True)
    opt_sdfs = [f"{mol_id}_opt.sdf" for mol_id in done_jobs_record.DFT_opt_freq if mol_id not in done_jobs_record.WFT_sp]
    for opt_sdf in opt_sdfs:
        mol_id = os.path.splitext(opt_sdf)[0].split("_")[0]
        os.makedirs(os.path.join(args.DLPNO_sp_folder, mol_id), exist_ok=True)
        shutil.copyfile(os.path.join(args.DFT_opt_freq_folder, mol_id, opt_sdf),
                        os.path.join(args.DLPNO_sp_folder, mol_id, mol_id + ".sdf"))
        charge = mol_id_to_charge_dict[mol_id]
        mult = mol_id_to_mult_dict[mol_id]
        os.chdir(os.path.join(args.DLPNO_sp_folder, mol_id))
        try:
            dlpno_sp_calc(mol_id, ORCA_PATH, charge, mult, args.DLPNO_sp_n_procs, args.DLPNO_sp_job_ram)
            done_jobs_record.WFT_sp.append(mol_id)
            done_jobs_record.save(project_dir, args.task_id)
            logger.info(f'DLPNO single point calculation for {mol_id} completed')
        except:
            logger.error(f'DLPNO single point calculation for {mol_id} failed.')
            logger.error(traceback.format_exc())
        os.chdir(project_dir)
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
    #         shutil.copyfile(os.path.join(args.semiempirical_opt_folder, opt_sdf),
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


