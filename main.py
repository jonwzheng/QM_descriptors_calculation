from argparse import ArgumentParser, Namespace
import os
import shutil

import pandas as pd

from lib import create_logger
from lib import csearch
from lib import xtb_optimization
from lib import dft_scf
from lib import dft_scf_dft_only

XTB_PATH = '/home/gridsan/jonzheng/.conda/envs/QM_descriptors/bin/'
G16_PATH = '/home/gridsan/groups/RMG/Software/gaussian/g16/'

parser = ArgumentParser()
parser.add_argument('--ismiles', type=str, required=False,
                    help='input smiles included in a .csv file')
parser.add_argument('--output', type=str, default='QM_descriptors.pickle',
                    help='output as a .pickle file')
# conformer searching
parser.add_argument('--MMFF_conf_folder', type=str, default='MMFF_conf',
                    help='folder for MMFF searched conformers')
parser.add_argument('--nconf', type=int, default=500,
                    help='number of MMFF conformers')
parser.add_argument('-max_conf_try', type=int, default=2000,
                    help='maximum attempt for conformer generating, '
                         'this is useful for molecules with many chiral centers.')
parser.add_argument('-rmspre', type=float, required=False,
                        help='rms threshold pre optimization')
parser.add_argument('--rmspost', type=float, required=False, default=0.4,
                    help='rms threshold post MMFF minimization')
parser.add_argument('--E_cutoff', type=float, required=False, default=10.0,
                    help='energy window for MMFF minimization')
parser.add_argument('--MMFF_threads', type=int, required=False, default=40,
                    help='number of process for the MMFF conformer searching')
parser.add_argument('--timeout', required=False, default=600,
                    help='time window for each MMFF conformer searching sub process')
# xtb optimization
parser.add_argument('--xtb_folder', type=str, default='XTB_opt',
                    help='folder for XTB optimization')

# DFT calculation
parser.add_argument('--DFT_folder', type=str, default='DFT',
                    help='folder for DFT calculation')
parser.add_argument('--DFT_theory', type=str, default='b3lyp/def2svp',
                    help='level of theory for the DFT calculation')
parser.add_argument('--DFT_n_procs', type=int, default=20,
                    help='number of process for DFT calculations')
parser.add_argument('--jobtype', type=str, default='neutral',
                    help='type of job (neutral, plus1, minus1). Default is neutral')
parser.add_argument('--DFT_only', type=bool, default=False,
                    help='Whether to only run a DFT calc using pre-specified chkpoint files')

args = parser.parse_args()


qm_descriptors = []
if args.DFT_only == False:
    name = os.path.splitext(args.ismiles)[0]
    logger = create_logger(name=name)

    df = pd.read_csv(args.ismiles, index_col=0)
    # conformer searching

    logger.info('starting MMFF conformer searching')
    supp = (x for x in df[['id', 'smiles']].values)
    conf_sdfs = csearch(supp, len(df), args, logger)

    # xtb optimization

    logger.info('starting GFN2-XTB structure optimization for the lowest MMFF conformer')
    if not os.path.isdir(args.xtb_folder):
        os.mkdir(args.xtb_folder)

    opt_sdfs = []
    for conf_sdf in conf_sdfs:
        try:
            shutil.copyfile(os.path.join(args.MMFF_conf_folder, conf_sdf),
                            os.path.join(args.xtb_folder, conf_sdf))
            opt_sdf = xtb_optimization(args.xtb_folder, conf_sdf, XTB_PATH, logger)
            opt_sdfs.append(opt_sdf)
        except Exception as e:
            logger.error('XTB optimization for {} failed: {}'.format(os.path.splitext(conf_sdf)[0], e))
    
    # G16 DFT calculation
    if not os.path.isdir(args.DFT_folder):
        os.mkdir(args.DFT_folder)

    for opt_sdf in opt_sdfs:
        try:
            shutil.copyfile(os.path.join(args.xtb_folder, opt_sdf),
                            os.path.join(args.DFT_folder, opt_sdf))
            qm_descriptor = dft_scf(args.DFT_folder, opt_sdf, G16_PATH, args.DFT_theory, args.DFT_n_procs,
                                    logger, args.jobtype)
            qm_descriptors.append(qm_descriptor)
        except Exception as e:
            logger.error('Gaussian optimization for {} failed: {}'.format(os.path.splitext(opt_sdf)[0], e))
else:
    logger = create_logger(name='dft')
    opt_chks = []
    # G16 DFT calculation
#    if not os.path.isdir(args.DFT_folder):
#        os.mkdir(args.DFT_folder)

    for file in os.listdir(os.getcwd()):
        if file.endswith(".chk"):
            opt_chks.append(file)

    for opt_chk in opt_chks:
        try:
            qm_descriptor = dft_scf_dft_only(args.DFT_folder, opt_chk, G16_PATH, args.DFT_theory, 
                                args.DFT_n_procs, logger, args.jobtype)
            qm_descriptors.append(qm_descriptor)
        except Exception as e:
            logger.error('Gaussian optimization for {} failed: {}'.format(os.path.splitext(opt_chk)[0], e))

qm_descriptors = pd.DataFrame(qm_descriptors)
qm_descriptors.to_pickle(args.output)
qm_descriptors.to_csv("output.csv")

