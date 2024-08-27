#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 15 12:48:51 2022

@author: ikalvet
"""
import argparse
import pyrosetta as pyr
import pyrosetta.rosetta
import pyrosetta.distributed.io
import sys, os
import itertools
import time
import numpy as np
import pandas as pd
import multiprocessing
import queue
import threading
import random
import scipy.spatial
script_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(script_dir)
sys.path.append(script_dir+'/utils/')
import protocol
import utils
import dunbrack_rotlib
import align_pdbs




def process_rotamer_set_queue(q, prefix, bad_rotamers, rotamers, cst_io, motifs):
    while True:
        _s = q.get()
        if _s is None:
            return

        i = _s[0]
        ids = _s[1]
        # Grabbing a combination of inverse rotamers based on the provided
        # per-cst inverse rotamer ids.
        c = [rotamers[n][i] for n, i in enumerate(ids)]

        if any([rot_id in bad_rotamers[j] for j, rot_id in enumerate(ids)]):
            # print(f"Bad rotamer in set {i}")
            continue

        # TODO: implement symmetry here
        # Take the list "c" and apply some symmetric transform to the residues there
        # Then the rest of the code should take care of it appropriately

        pose = pyrosetta.rosetta.core.pose.Pose()
        bad_rotamer = False
        catres_resnos = {n: 0 for n,r in enumerate(c) if not isinstance(r, pyrosetta.rosetta.core.pose.Pose) and r.is_ligand()}
        ligands = [r for r in c if not isinstance(r, pyrosetta.rosetta.core.pose.Pose) and r.is_ligand()]
        for j, res in enumerate(c):
            if args.debug:
                if isinstance(res, pyrosetta.rosetta.core.conformation.Residue):
                    print(i, j, res.name())
                elif isinstance(res, pyrosetta.rosetta.core.pose.Pose):
                    print(i, j, res.pdb_info().name())

            if not isinstance(res, pyrosetta.rosetta.core.pose.Pose) and res.is_ligand():  # ligand
                continue

            # If we have already seen that it's a bad rotamer then let's just skip it
            if ids[j] in bad_rotamers[j]:
                if args.debug: print(f"{j}, previously seen as a bad rotamer")
                bad_rotamer = True
                break

            if isinstance(res, pyrosetta.rosetta.core.conformation.Residue):
                _res_pose = pyrosetta.rosetta.core.pose.Pose()
                _res_pose.append_residue_by_jump(res, 0)
                if res.is_protein():
                    _res_pose = protocol.extend_SS(pose=_res_pose, ref_seqpos=1,
                                          secstruct=args.secstruct_per_cst[j], AAA=AAA,
                                          nres_Nterm=args.N_len_per_cst[j],
                                          nres_Cterm=args.C_len_per_cst[j])
                    _res_pose.fold_tree().clear()
                    _res_pose.fold_tree().add_edge(1, _res_pose.size(), -1)  # This will avoid FoldTree reordering error showing up
                catres_resno = args.N_len_per_cst[j]+1

            elif isinstance(res, pyrosetta.rosetta.core.pose.Pose):
                _res_pose = res.clone()
                catres_resno = motifs[j]["resno"]


            # Adding ligand to the extended chain and checking for clashes
            for ligand in ligands:
                # _res_pose.append_residue_by_jump(ligand, 1)  # this doesn't turn ligand into new chain
                _res_pose.append_residue_by_jump(ligand, catres_resno,
                                                 jump_anchor_atom=_res_pose.residue(catres_resno).atom_name(_res_pose.residue(catres_resno).nbr_atom()),
                                                 jump_root_atom=ligand.atom_name(ligand.nbr_atom()),
                                                 start_new_chain=True)

            if protocol.check_clash(_res_pose, catres_resnos=[catres_resno]+[r.seqpos() for r in _res_pose.residues if r.is_ligand()], tip_atom=args.tip_atom, debug=args.debug) is True:
                if args.debug: print(f"{j}, clash after extension")
                # Only adding the residude object to the bad residues
                # The motif pose will never be dumped
                if isinstance(res, pyrosetta.rosetta.core.conformation.Residue):
                    bad_rotamers[j].append(ids[j])
                    # print(j, ids[j], list(bad_rotamers[j]))
                elif isinstance(res, pyrosetta.rosetta.core.pose.Pose):
                    if args.debug: print("MOTIF POSE SEEMS TO GIVE CLASH!!!! PLEASE INVESTIGATE!!!")
                bad_rotamer = True

                # Giving up if all rotamers are bad
                if len(bad_rotamers[j]) == len(rotamers[j]):
                    print(f"All rotamers for CST {j} are bad...")
                break

            if isinstance(res, pyrosetta.rosetta.core.conformation.Residue):
                catres_resnos[j] = pose.size() + args.N_len_per_cst[j]+1
            else:
                catres_resnos[j] = motifs[j]["resno"]

            pyrosetta.rosetta.core.pose.append_subpose_to_pose(pose, _res_pose, 1, _res_pose.size()-len(ligands), new_chain=True)

        # Finished individual evaluation of residues
        # Now putting the whole thing together
        if bad_rotamer is True:
            if args.debug: print(f"{j}, bad rotamer")
            continue

        # Adding ligand as the last residue
        for _n,res in enumerate(c):
            if isinstance(res, pyrosetta.rosetta.core.pose.Pose):
                continue
            if res.is_ligand():
                lig_pose = pyrosetta.rosetta.core.pose.Pose()
                lig_pose.append_residue_by_jump(res, 0)
                pyrosetta.rosetta.core.pose.append_subpose_to_pose(pose, lig_pose, 1, 1, new_chain=True)
                catres_resnos[_n] = pose.size()

        # Checking for clashes
        # Ignoring clashes between catalytic residues and the ligand
        ignore_clash_respairs = []
        for j in catres_resnos:
            if isinstance(c[j], pyrosetta.rosetta.core.conformation.Residue):
                assert pose.residue(catres_resnos[j]).name3() == c[j].name3(), f"cst {j}: resno {catres_resnos[j]}, {c[j].name3()} != {pose.residue(catres_resnos[j]).name3()}"
            if j == 0:
                continue
            if args.debug: print(f"clashcheck exclude cst atoms, cst {j}, resno {catres_resnos[j]}, name {pose.residue(catres_resnos[j]).name()}")
            ignore_clash_respairs.append((catres_resnos[0], catres_resnos[j]))

        clash = protocol.check_clash(pose, catres_resnos=catres_resnos.values(), ignore_respairs=ignore_clash_respairs, tip_atom=args.tip_atom, debug=args.debug)
        if clash is True:
            if args.debug: print(f"{j}, clash in the final assembly")
            continue
        if args.debug: print(j, pose.sequence())
        
        # TODO: Need to implement checking whether the pose actually respects the CST's
        # This is an issue when the ligand has any chi sampling enabled, and another residue is matched downstream of that.
        # Some combinations of rotamers are not meant to work together
        ## I think this is now managed in the REMARK 666 generation stage

        pose_name = args.prefix
        for res in c:
            if isinstance(res, pyrosetta.rosetta.core.conformation.Residue):
                if res.is_protein():
                    pose_name += res.name1() + "_"
                else:
                    pose_name += res.name3() + "_"
            elif isinstance(res, pyrosetta.rosetta.core.pose.Pose):
                pose_name += os.path.basename(res.pdb_info().name()).replace(".pdb", "") + "_"
        pose_name += f"{prefix}_{i}{args.suffix}.pdb"
        if os.path.exists(pose_name):
            print(f"Found existing file with name {pose_name}")
            pose_name.replace(".pdb", "a.pdb")


        remarks = protocol.create_remark_lines(pose, catres_resnos, cst_io)

        if len(remarks) != len(catres_resnos) - 1:
            if args.debug: print(f"{i}: Could not build all REMARK 666 lines")
            continue

        print(f"Found good rotamer: {pose_name.replace('.pdb', '')}")

        pdbstr = pyrosetta.distributed.io.to_pdbstring(pose).split("\n")

        pdbstr_new = []
        for l in pdbstr:
            pdbstr_new.append(l)
            if "HEADER" in l:
                for rmrk in remarks:
                    pdbstr_new.append(rmrk)
        with open(pose_name, "w") as file:
            file.write("\n".join(pdbstr_new))



def parallelize_mp(iterables, rotset, prefix, cst_io, motifs):
    print(f"{len(iterables)} configurations to process")
    the_queue = multiprocessing.Queue()  # Queue stores the iterables

    start = time.time()
    manager = multiprocessing.Manager() 
    bad_rotamers = manager.dict()  # Need a special dictionary to store outputs from multiple processes
    results_found = manager.dict()

    for i, c in enumerate(iterables):
        the_queue.put((i, c))

    for j in range(len(c)):
        bad_rotamers[j] = manager.list()

    print(f"Starting to generate inverse rotamer assemblies using {args.nproc} parallel processes.")
    pool = multiprocessing.Pool(processes=args.nproc,
                                initializer=process_rotamer_set_queue,
                                initargs=(the_queue, prefix, bad_rotamers, rotset, cst_io, motifs, ))

    # None to end each process
    for _i in range(args.nproc):
        the_queue.put(None)

    # Closing the queue and the pool
    the_queue.close()
    the_queue.join_thread()
    pool.close()
    pool.join()
    
    for j in bad_rotamers:
        print(j, list(bad_rotamers[j]))

    end = time.time()
    print(f"Processing all the rotamers in set {prefix} took {(end - start):.2f} seconds")




def main(args):
    if args.suffix != "":
        args.suffix = f"_{args.suffix}"

    if args.prefix != "":
        args.prefix = f"{args.prefix}"

    assert os.path.exists(args.cstfile)
    extra_res_fa = ""
    if args.params is not None:
        params = [p for p in args.params if ".params" in p]
        extra_res_fa = "-extra_res_fa " + ' '.join(params)

    """
    Setting up PyRosetta
    """
    
    # pyr.init(f"{extra_res_fa} -run:preserve_header -output_virtual true")
    pyr.init(f"{extra_res_fa} -run:preserve_header")
    
    # Loading the backbone-dependent Dunbrack rotamer library into a dataframe
    dunbrack_database = os.path.dirname(pyr.__file__) + "/database/rotamer/bbdep02.May.sortlib-correct.12.2010"
    rotlib = dunbrack_rotlib.load_rotamer_df(dunbrack_database)


    global AAA  # making it global so that functions downstream can see it
    AAA = pyr.pose_from_sequence("AAA")


    ###### CST PARSING ########
    # Parsing the CST file
    addcst_mover = pyrosetta.rosetta.protocols.enzdes.AddOrRemoveMatchCsts()
    chem_manager = pyrosetta.rosetta.core.chemical.ChemicalManager.get_instance()
    residue_type_set = chem_manager.residue_type_set("fa_standard")
    cst_io = pyrosetta.rosetta.protocols.toolbox.match_enzdes_util.EnzConstraintIO(residue_type_set)
    cst_io.read_enzyme_cstfile(args.cstfile)


    # Figuring out which residue atoms are used for each cst
    # Using the MCFI (MatcherConstraintFileInfo) object for that
    # cst_atoms will be a dict where each cst_block contains a list of variable CST's? and then a list of residue types
    cst_atoms = protocol.get_cst_atoms(cst_io)


    # Storing information about which residues are matched for each CST block
    restypes = {}
    for n in range(1, cst_io.mcfi_lists_size()+1):
        restypes[n] = []
        for restype in cst_io.mcfi_list(n).upstream_restypes():
            restypes[n].append(restype.name3())


    ### PROCESS ARGUMENTS A BIT FURTHER ###
    args = protocol.parse_arguments(args, restypes)


    #### PARSING HIS TAUTOMER RESTRICTIONS #####
    keep_his_tautomer_per_cst = None
    if args.keep_his_tautomer is not None:
        keep_his_tautomer_per_cst = {int(x.split(":")[0]): x.split(":")[1] for x in args.keep_his_tautomer.split(",")}
        assert all([val in ["HIS", "HIS_D"] for key, val in keep_his_tautomer_per_cst.items()]), "Invalid input for --keep_his_tautomer"


    ### ROTAMER SUBSAMPLING ####
    chi_subsampling_levels = protocol.parse_rotamer_subsampling(args, cst_atoms)


    ### Putting together a dictionary listing good rotamers for each residue in each CST
    restype_good_rotamers = {}
    for n in restypes:
        restype_good_rotamers[n] = {}
        for restyp in restypes[n]:
            if restyp not in utils.N_chis.keys():
                continue
            if restyp not in restype_good_rotamers.keys():
                use_only_best_rotamer = False
                if n in args.use_best_rotamer_cstids:
                    use_only_best_rotamer = True
                restype_good_rotamers[n][restyp] = dunbrack_rotlib.find_good_rotamers(rotlib, restyp, args.dunbrack_prob_per_cst[n],
                                                                                      args.secstruct_per_cst[n],
                                                                                      keep_only_best=use_only_best_rotamer)


    ### PARSING EXTERNAL MOTIFS ####
    # TODO: make external motifs usable with other CST id's, not just the 1st one
    motifs = None
    if args.motif_for_cst is not None:
        motifs = protocol.parse_motif_input(args.motif_for_cst, cst_atoms, restypes)

    
    
    ### GETTING INVERSE ROTAMERS ####
    ### This is where half of the work gets done ###
    invrot_tree = pyrosetta.rosetta.protocols.toolbox.match_enzdes_util.TheozymeInvrotTree(cst_io)
    invrot_tree.generate_targets_and_inverse_rotamers()
    all_inverse_rotamers_per_cst = invrot_tree.collect_all_inverse_rotamers()
    
    
    ## There is a way to get inverse rotamers from cst_io
    ## need to investigate this, because this allows keeping the sub-cst information
    """
    target_ats = pyrosetta.rosetta.utility.vector1_unsigned_long()
    invrot_ats = pyrosetta.rosetta.utility.vector1_unsigned_long()
    
    _mcfi.inverse_rotamers_against_residue(target_conf=lig, invrot_restype=_mcfi.allowed_restypes(_mcfi.upstream_res())[1],
                                           target_ats=target_ats, invrot_ats=invrot_ats, flip_exgs_upstream_downstream_samples=False, backbone_interaction=False)
    """
    

    time.sleep(1)
    
    print(f"{len(all_inverse_rotamers_per_cst)} rotamer sets to process")
    for xx, rotset in enumerate(all_inverse_rotamers_per_cst):
        print(f"Non-redundant rotamer set {xx+1}")
        for cst_block, invrots in enumerate(rotset.invrots()):
            print(f"CST {cst_block}: {len(invrots)} inverse rotamers.")

        # Listify the inverse rotamer dataset
        rotset_sub = [[invrot for invrot in invrots] for invrots in rotset.invrots()]

        # Pruning all other inverse rotamers based on proton-chis.
        # Removing duplicate rotamers where the only difference is in the value of the proton_chi
        for rotset_id in range(len(rotset_sub)):
            if isinstance(rotset_sub[rotset_id][0], pyrosetta.rosetta.core.pose.Pose) or rotset_sub[rotset_id][0].is_ligand():
                continue
            _n_before = len(rotset_sub[rotset_id])
            rotset_sub[rotset_id] = protocol.prune_residue_rotamers(rotset_sub[rotset_id])
            if len(rotset_sub[rotset_id]) != _n_before:
                print(f"CST {rotset_id}: {len(rotset_sub[rotset_id])} inverse rotamers after pruning for proton-chi")


        # Loading any external motifs, if provided and aligning them to the appropriate CST atoms
        if args.motif_for_cst is not None:
            for cstno in motifs:
                # TODO: implement for not-first CST's (or CST's with additional sampling from CST file),
                # Picking rotamers with unique subsampling defined in CST
                to_align_rotamers = protocol.find_unique_rotamers_for_motif([r if i==cstno else [] for i, r in enumerate(rotset_sub)], motifs)
                rotset_sub[cstno] = [align_pdbs.align_pose_to_residue(rotamer, motifs[cstno]["pose"],
                                                                     {"atoms1": motifs[cstno]["atoms"],
                                                                      "atoms2": [(motifs[cstno]["resno"], a) for a in motifs[cstno]["atoms"]]}) for rotamer in to_align_rotamers[cstno]]


        # Pruning inverse rotamers based on Dunbrack probabilites
        rotset_sub = protocol.preselect_inverse_rotamers(rotset_sub, restype_good_rotamers, keep_his_tautomer_per_cst)
        if rotset_sub is None:
            continue

        # Culling ligand rotamers based on RMSD cutoff
        if args.prune_ligand_rotamers != 0.0:
            for rotset_id in range(len(rotset_sub)):
                if isinstance(rotset_sub[rotset_id][0], pyrosetta.rosetta.core.pose.Pose):
                    continue
                if rotset_sub[rotset_id][0].is_ligand():
                    rotset_sub[rotset_id] = protocol.prune_ligand_rotamers(rotset_sub[rotset_id], args.prune_ligand_rotamers, args.nproc)

        # Performing rotamer subsampling (expanding CHI's)
        if any([any([y != 0 for y in x.values()]) for k, x in chi_subsampling_levels.items()]):
            rotset_sub = protocol.subsample_rotamers(rotset_sub, chi_subsampling_levels, restype_good_rotamers, cst_atoms)

        # Picking random rotamers if requested
        if args.frac_random_rotamers_per_cst is not None or args.max_random_rotamers_per_cst is not None:
            print("Picking a random subset of inverse rotamers")
            rotset_sub = protocol.pick_random_rotamers_set(rotset_sub, max_random_rotamers_per_cst=args.max_random_rotamers_per_cst,
                                                  frac_random_rotamers_per_cst=args.frac_random_rotamers_per_cst)
    
        for cst_block, invrots in enumerate(rotset_sub):
            print(f"CST {cst_block}: {len(invrots)} inverse rotamers after filtering.")
    
        rotset_ids = [[i for i, y in enumerate(x)] for x in rotset_sub]
        combs = itertools.product(*[x for x in rotset_ids])

        # Processing this subset of rotamers
        parallelize_mp(iterables=[c for c in combs], rotset=rotset_sub, prefix=xx+1, cst_io=cst_io, motifs=motifs)



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--cstfile", type=str, required=True, help="CST file used for matching. Keep sampling to minimum to avoid combinatorial explosion.")
    parser.add_argument("--params", nargs="+", required=False, help="params files used by ligands and residues")
    parser.add_argument("--keep_his_tautomer", type=str, help="Per cst, should a specific HIS tautomer (HIS or HIS_D) be used. Keeps only one the requested HIS tautomers. Format: 'cst_no:HIS/HIS_D,..'")
    parser.add_argument("--dunbrack_prob", type=float, default=0.85, help="Cumulative Dunbrack probability of used rotamers for any residue\n."
                                                                          "As used by the -packing:dunbrack_prob_... flag in Rosetta.")
    parser.add_argument("--dunbrack_prob_per_cst", type=float, nargs="+", help="Cumulative Dunbrack probability of used rotamers for each CST residue.")
    parser.add_argument("--N_len", type=int, default=4, help="Number of residues added to the stub N-term")
    parser.add_argument("--C_len", type=int, default=5, help="Number of residues added to the stub C-term")
    parser.add_argument("--N_len_per_cst", type=int, nargs="+", help="Number of residues added to the stub N-term, per CST")
    parser.add_argument("--C_len_per_cst", type=int, nargs="+", help="Number of residues added to the stub C-term, per CST")
    parser.add_argument("--prune_ligand_rotamers", type=float, default=0.0, help="Prunes the set of used ligand rotamers based on clashcheck, AND rmsd similarity cutoff.")
    parser.add_argument("--max_random_rotamers", type=int, help="Number of random rotamers picked for each residue for the sampling. Reasonable number would be below 20 for quick sampling.")
    parser.add_argument("--max_random_rotamers_per_cst", nargs="+", type=int, help="Number of random rotamers picked for each CST block for the sampling. First value is for the ligand.")
    parser.add_argument("--frac_random_rotamers", type=float, help="Fraction of rotamers that are randomly picked for each residue for the sampling.")
    parser.add_argument("--frac_random_rotamers_per_cst", nargs="+", type=float, help="Fraction of rotamers that are randomly picked for each CST block for the sampling. First value is for the ligand.")
    parser.add_argument("--secstruct", type=str, default="H", choices=["E", "H"], help="What secondary structure stub should be generated for each residue.")
    parser.add_argument("--secstruct_per_cst", nargs="+", type=str, help="Per CST, what secondary structure stub should be generated for each residue.")
    parser.add_argument("--motif_for_cst", type=str, nargs="+", help="Per CST, an external motif that should be used, instead of inverse rotamers. Only works for the first CST right now. Format: cst_no:resno_in_motif:filepath ...")
    parser.add_argument("--use_best_rotamer_cstids", nargs="+", type=int, default=[], help="CST ID's that should only use the best rotamer from each secondary structure bin. Numbering starts from 1.")
    parser.add_argument("--extra_chi", type=str, help="Enables extra CHI sampling on a given level for all CST's. Input format: chino:level,chino2:level2")
    parser.add_argument("--extra_chi_per_cst", nargs="+", help=f"Enables extra CHI sampling on a given level for specific CST's. Input format: CSTNO1-chino:level,chino2:level2 CSTNO2-chino:level,chino2:level2\nSampling levels:\n{protocol.calculate_samplings.__doc__}")
    parser.add_argument("--suffix", type=str, default= "", help="Suffix to be added to the end of output files")
    parser.add_argument("--prefix", type=str, default= "", help="Prefix to be added to the beginning of output files")
    parser.add_argument("--tip_atom", action="store_true", default=False, help="Inverse rotamers will be pre-selected based on whether the tip atoms are placed geometrically differently. Rotamer diversity is ignored.")
    parser.add_argument("--nproc", type=int, help="Number of CPU cores used.")
    parser.add_argument("--debug", action="store_true", default=False, help="Debug mode. Will print out more output at each step. Will run in single-core mode.")

    args = parser.parse_args()

    if "SLURM_CPUS_ON_NODE" in os.environ:
        args.nproc = int(os.environ["SLURM_CPUS_ON_NODE"])
    if args.nproc is None:
        args.nproc = os.cpu_count()
    if args.debug is True:
        args.nproc = 1

    main(args)

