import logging
import os
import json
from rdkit import Chem

def create_logger(name: str, task_id: int) -> logging.Logger:
    """
    Creates a logger with a stream handler and two file handlers.

    The stream handler prints to the screen depending on the value of `quiet`.
    One file handler (verbose.log) saves all logs, the other (quiet.log) only saves important info.

    :param save_dir: The directory in which to save the logs.
    :return: The logger.
    """
    logging.basicConfig(
            filemode='w+',
            level=logging.INFO)
    logger = logging.getLogger(name)
    logger.propagate = False
    file_name = f'{name}_{task_id}.log'
    try:
        os.remove(file_name)
    except:
        pass
    fh = logging.FileHandler(filename=file_name)
    fh.setLevel(logging.DEBUG)
    logger.addHandler(fh)

    return logger

def write_mol_to_sdf(mol, path, confIds=[0], confEns=None):
    if isinstance(confIds, int):
        confIds = [confIds]
    if isinstance(confEns, int):
        confEns = [confEns]
    writer = Chem.SDWriter(path)
    if confEns:
        for confId, confEn in zip(confIds, confEns):
            mol.SetProp('ConfId', str(confId))
            mol.SetProp('ConfEnergies', str(confEn) + ' kcal/mol')
            writer.write(mol, confId=confId)
    else:
        for confId in confIds:
            writer.write(mol, confId=confId)
    writer.close()

def load_sdf(path, removeHs=False, sanitize=False):
    return Chem.SDMolSupplier(path, removeHs=removeHs, sanitize=sanitize)

class DoneJobsRecord(object):
    """
    class to record completed jobs
    """
    def __init__(self):
        self.FF_conf = []
        self.XTB_opt_freq = []
        self.DFT_opt_freq = []
        self.COSMO = {}
        self.WFT_sp = []
        self.QM_desp = []
    
    def save(self, project_dir, task_id):
        with open(os.path.join(project_dir, f"done_jobs_record_{task_id}.json"), "w+") as fh:
            json.dump(vars(self), fh)

    def load(self, project_dir, task_id):
        with open(os.path.join(project_dir,f"done_jobs_record_{task_id}.json"), "r") as fh:
            content = json.load(fh)
        for job, molids in content.items():
            setattr(self, job, molids)

done_jobs_record = DoneJobsRecord()