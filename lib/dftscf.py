from rdkit import Chem
import os
import subprocess
from .file_parser import mol2xyz, xyz2com
from .grab_QM_descriptors import read_log


def dft_scf(folder, sdf, g16_path, level_of_theory, n_procs, logger, jobtype):
    basename = os.path.basename(sdf)

    parent_folder = os.getcwd()
    os.chdir(folder)

    try:
        file_name = os.path.splitext(basename)[0]

        xyz = mol2xyz(Chem.SDMolSupplier(sdf, removeHs=False, sanitize=False)[0])

        pwd = os.getcwd()

        g16_command = os.path.join(g16_path, 'g16')
        QM_descriptors = {}

        if not os.path.isdir(jobtype):
            os.mkdir(jobtype)

        if jobtype == 'neutral':
            charge = 0
            mult = 1
        elif jobtype == 'plus1':
            charge = 1
            mult = 1
        elif jobtype == 'minus1':
            charge = -1
            mult = 1

        head = '%chk={}.chk\n%nprocshared={}\n# {} scf=(maxcycle=512, xqc) opt freq iop(2/9=2000)' \
                '\n'.format(file_name, n_procs, level_of_theory)

        os.chdir(jobtype)
        comfile = file_name + '.gjf'
        xyz2com(xyz, head=head, comfile=comfile, charge=charge, mult=mult, footer='$NBO BNDIDX $END\n')

        logfile = file_name + '.log'
        outfile = file_name + '.out'
        with open(outfile, 'w') as out:
            subprocess.run('{} < {} >> {}'.format(g16_command, comfile, logfile), shell=True, stdout=out, stderr=out)
            QM_descriptors[jobtype] = read_log(logfile)
        os.chdir(pwd)

        QM_descriptor_return = QM_descriptors[jobtype]

        os.remove(sdf)
    finally:
        os.chdir(parent_folder)

    return QM_descriptor_return


def dft_scf_dft_only(folder, chk, g16_path, level_of_theory, n_procs, logger, jobtype):
    basename = os.path.basename(chk)

    parent_folder = os.getcwd()


    try:
        file_name = os.path.splitext(basename)[0]


        pwd = os.getcwd()

        g16_command = os.path.join(g16_path, 'g16')
        QM_descriptors = {}

        if not os.path.isdir(jobtype):
            os.mkdir(jobtype)

        if jobtype == 'neutral':
            charge = 0
            mult = 1
        elif jobtype == 'plus1':
            charge = 1
            mult = 1
        elif jobtype == 'minus1':
            charge = -1
            mult = 1

        head = '%chk={}.chk\n%nprocshared={}\n# {} scf=(maxcycle=512, xqc) sp Guess=Read Geom=AllCheckpoint iop(2/9=2000)' \
                '\n'.format(file_name, n_procs, level_of_theory)

        comfile = file_name + '.gjf'
        with open(comfile, 'w') as com:
            com.write(head)

        logfile = file_name + '.log'
        outfile = file_name + '.out'
        with open(outfile, 'w') as out:
            subprocess.run('{} < {} >> {}'.format(g16_command, comfile, logfile), shell=True, stdout=out, stderr=out)
            QM_descriptors[jobtype] = read_log(logfile)

        QM_descriptor_return = QM_descriptors[jobtype]

    finally:
        os.chdir(parent_folder)

    return QM_descriptor_return
