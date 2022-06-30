#!/usr/bin/python
from __future__ import print_function, absolute_import
import shutil
import subprocess

from multiprocessing import Process, Manager

from rdkit import Chem
from rdkit.Chem import AllChem
from concurrent import futures
from .log_parser import XtbLog
from .file_parser import write_mol_to_sdf, load_sdf
import os
import traceback

# algorithm to generate nc conformations
def _genConf(s, max_n_conf, max_try, rms, E_cutoff_fraction, rmspost, return_dict, name, conf_search_FF, FF_conf_folder, XTB_path):
    m = Chem.MolFromSmiles(s)
    if not m:
        return
    try:
        AllChem.EmbedMolecule(m, AllChem.ETKDG())
        m = Chem.AddHs(m, addCoords=True)
    except:
        return

    nr = int(AllChem.CalcNumRotatableBonds(m))

    tnr = 3**nr
    nc = tnr if tnr < max_n_conf else max_n_conf

    if not rms:
        rms = -1
    ids = AllChem.EmbedMultipleConfs(m, numConfs=nc, maxAttempts=max_try, pruneRmsThresh=rms,
                                   randomSeed=1, useExpTorsionAnglePrefs=True, useBasicKnowledge=True)

    if len(ids)== 0:
        ids = m.AddConformer(m.GetConformer(), assignID=True)

    diz = []
    pwd = os.getcwd()

    try:
        for id in ids:
            if conf_search_FF == "MMFF94s":
                prop = AllChem.MMFFGetMoleculeProperties(m, mmffVariant="MMFF94s")
                ff = AllChem.MMFFGetMoleculeForceField(m, prop, confId=id)
                ff.Minimize()
                en = float(ff.CalcEnergy())
            elif conf_search_FF == "GFNFF":
                scratch_dir = os.path.join(FF_conf_folder, f'{name}_{id}')
                os.makedirs(scratch_dir, exist_ok=True)
                os.chdir(scratch_dir)

                input_file_name = f'{name}_{id}.sdf'
                write_mol_to_sdf(m, input_file_name, id)

                xtb_command = os.path.join(XTB_path, 'xtb')
                output_file_name = f'{name}_{id}.log'
                with open(output_file_name, 'w') as out:
                    subprocess.call([xtb_command, '--gfnff', input_file_name, '--opt'],
                                    stdout=out, stderr=out)

                log = XtbLog(output_file_name)
                en = float(log.E)*627.5 #convert to kcal/mol
                opt_conf = load_sdf("xtbopt.sdf")[0].GetConformer()
                conf = m.GetConformer(id)
                for i in range(m.GetNumAtoms()):
                    pt = opt_conf.GetAtomPosition(i)
                    conf.SetAtomPosition(i, (pt.x, pt.y, pt.z))
                os.chdir(pwd)
                #shutil.rmtree(scratch_dir)
            econf = (en, id)
            diz.append(econf)
    except Exception as e:
        print(traceback.format_exc())
        return_dict['return'] = (None, None, None)
        return
    
    if E_cutoff_fraction != "Y":
        n, diz2 = energy_filter(m, diz, E_cutoff_fraction)
    else:
        n = m
        diz2 = diz

    if rmspost is not None and n.GetNumConformers() > 1:
        o, diz3 = postrmsd(n, diz2, rmspost)
    else:
        o = n
        diz3 = diz2
    return_dict['return'] = (o, diz3, nr)


# wrap the genConf in process so that the genConf can be stopped
class genConf:
    def __init__(self, m, args):
        chembl_id, SMILES = m
        self.s = SMILES
        self.name = chembl_id
        self.max_n_conf = args.max_n_conf
        self.max_nc_try = args.max_conf_try
        self.rms = args.rmspre
        self.E_cutoff_fraction = args.E_cutoff_fraction
        self.rmspost = args.rmspost
        self.timeout = args.timeout
        self.conf_search_FF = args.conf_search_FF
        self.FF_conf_folder = args.FF_conf_folder
        self.XTB_path = args.XTB_path
        
    def __call__(self):
        self.return_dict = Manager().dict()
        self.process = Process(target=_genConf, args=(self.s, self.max_n_conf, self.max_nc_try, self.rms, self.E_cutoff_fraction,
                                                      self.rmspost, self.return_dict, self.name, self.conf_search_FF, self.FF_conf_folder, self.XTB_path))

        self.process.start()
        self.process.join(self.timeout)
        if 'return' in self.return_dict:
            return self.return_dict['return']
        else:
            self.terminate()
            return (None, None, None)

    def terminate(self):
        self.process.terminate()


# filter conformers based on relative energy
def energy_filter(m, diz, E_cutoff_fraction):
    diz.sort()
    mini = float(diz[0][0])
    sup = mini + abs(mini) * E_cutoff_fraction
    n = Chem.Mol(m)
    n.RemoveAllConformers()
    n.AddConformer(m.GetConformer(int(diz[0][1])))
    nid = []
    ener = []
    nid.append(int(diz[0][1]))
    ener.append(float(diz[0][0])-mini)
    del diz[0]
    for x,y in diz:
        if x <= sup:
            n.AddConformer(m.GetConformer(int(y)))
            nid.append(int(y))
            ener.append(float(x-mini))
        else:
            break
    diz2 = list(zip(ener, nid))
    return n, diz2


# filter conformers based on geometric RMS
def postrmsd(n, diz2, rmspost):
    diz2.sort(key=lambda x: x[0])
    o = Chem.Mol(n)
    confidlist = [diz2[0][1]]
    enval = [diz2[0][0]]
    nh = Chem.RemoveHs(n)
    del diz2[0]
    for z,w in diz2:
        confid = int(w)
        p=0
        for conf2id in confidlist:
            rmsd = AllChem.GetBestRMS(nh, nh, prbId=confid, refId=conf2id)
            if rmsd < rmspost:
                p=p+1
                break
        if p == 0:
            confidlist.append(int(confid))
            enval.append(float(z))
    diz3 = list(zip(enval, confidlist))
    return o, diz3


# conformational search / handles parallel threads if more than one structure is defined
def csearch(supp, total, args, logger, done_jobs_record, project_dir):
    os.makedirs(args.FF_conf_folder, exist_ok=True)
    try:
        with futures.ProcessPoolExecutor(max_workers=args.FF_threads) as executor:
            n_tasks = args.FF_threads if args.FF_threads < total else total
            tasks = [genConf(next(supp), args) for m in range(n_tasks)]
            running_pool = {task.name: executor.submit(task) for task in tasks}

            while True:
                if len(running_pool) == 0:
                    break

                for mol_id in list(running_pool):
                    future = running_pool[mol_id]
                    if future.done():
                        mol, ids, nr = future.result(timeout=0)
                        if mol:
                            lowest_en, lowest_id = ids[0]
                            mol.SetProp('_Name', mol_id)
                            os.makedirs(os.path.join(args.FF_conf_folder, mol_id), exist_ok=True)
                            write_mol_to_sdf(mol, os.path.join(args.FF_conf_folder, mol_id, '{}.sdf'.format(mol_id)), lowest_id, lowest_en)

                            conformers_found = len(ids)
                            ids_to_save = [id for (en, id) in ids[:args.n_lowest_E_confs_to_save]]
                            ens_to_save = [en for (en, id) in ids[:args.n_lowest_E_confs_to_save]]
                            logger.info('conformer searching for {} completed: '
                                        '{} conformers found, save the lowest {}'.format(mol_id, conformers_found, len(ids_to_save)))
                            write_mol_to_sdf(mol, os.path.join(args.FF_conf_folder, mol_id, '{}_confs.sdf'.format(mol_id)), ids_to_save, ens_to_save)
                            done_jobs_record.FF_conf.append(mol_id)
                            done_jobs_record.save(project_dir, args.task_id)
                        else:
                            logger.info('conformer searching for {} failed.'.format(mol_id))
                            pass

                        # add new task
                        del(running_pool[mol_id])
                        
                        try:
                            task = genConf(next(supp), args)
                        except StopIteration:
                            # reach end of the supp
                            pass
                        else:
                            running_pool[task.name] = executor.submit(task)
    except StopIteration:
        # supp is empty
        pass
    logger.info('FF conformer searching finished')
    return done_jobs_record
