#!/usr/bin/env python
# coding: utf-8

import os
from rdmc.mol import RDKitMol

def ff_conf_parser(mol_id, mol_smi, mol_confs_sdf=None):

    failed_job = dict()
    valid_job = dict()

    if mol_confs_sdf is None:
        ids = str(int(int(mol_id.split("id")[1])/1000)) 
        mol_confs_sdf = os.path.join("output", "FF_conf", "outputs", f"outputs_{ids}", f"{mol_id}_confs.sdf")

    if os.path.isfile(mol_confs_sdf):

        failed_job[mol_id] = dict()
        valid_job[mol_id] = dict()

        pre_adj = RDKitMol.FromSmiles(mol_smi).GetAdjacencyMatrix()
        
        mols = RDKitMol.FromFile(mol_confs_sdf)
        for conf_id, mol in enumerate(mols):
            post_adj = mol.GetAdjacencyMatrix()
            try:
                (pre_adj == post_adj).all()
            except:
                print(mol_confs_sdf)
                break
            
            if (pre_adj == post_adj).all():
                valid_job[mol_id][conf_id] = {}
                xyz = mol.ToXYZ()
                en = mol.GetProp("ConfEnergies")
                valid_job[mol_id][conf_id]["ff_xyz"] = xyz
                valid_job[mol_id][conf_id]["ff_energy"] = en
            else:
                failed_job[mol_id][conf_id] = 'adjacency matrix'
        
        if not valid_job[mol_id]:
            del valid_job[mol_id]
            failed_job[mol_id]['reason'] = 'all confs failed'
        else:
            if not failed_job[mol_id]:
                del failed_job[mol_id]
    else:
        failed_job[mol_id] = dict()
        failed_job[mol_id]['reason'] = 'file not found'

    return failed_job, valid_job