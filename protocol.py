#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 25 23:14:34 2024

@author: indrek
"""
import pyrosetta as pyr
import pyrosetta.rosetta
import os, sys
import random
import numpy as np
import itertools
import multiprocessing
import time
import scipy.spatial

script_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(script_dir)
sys.path.append(script_dir+'/utils/')
import utils
import dunbrack_rotlib
import align_pdbs


"""
PARSING FUNCTIONS
"""
def parse_arguments(args, restypes):
    # Limiting Dunbrack library as requested. 
    if args.dunbrack_prob_per_cst is None:
        args.dunbrack_prob_per_cst = [None]+[args.dunbrack_prob for r in restypes]
    else:
        assert all([isinstance(x, float) for x in args.dunbrack_prob_per_cst])
        args.dunbrack_prob_per_cst = [None]+ args.dunbrack_prob_per_cst


    ######### IF REQUESTED... ############
    ### RANDOM ROTAMER SELECTION SETUP ###
    if args.max_random_rotamers_per_cst is not None:
        assert all([isinstance(x, int) for x in args.max_random_rotamers_per_cst])
        assert len(args.max_random_rotamers_per_cst) == len(restypes)+1, "Invalid number of per-cst max_random_rotamers_per_cst"
    
    if args.frac_random_rotamers_per_cst is not None:
        assert all([isinstance(x, float) for x in args.frac_random_rotamers_per_cst])
        assert len(args.frac_random_rotamers_per_cst) == len(restypes)+1, "Invalid number of per-cst frac_random_rotamers_per_cst"
    
    if args.max_random_rotamers is not None:
        args.max_random_rotamers_per_cst = [args.max_random_rotamers]+[args.max_random_rotamers for r in restypes]
    
    if args.frac_random_rotamers is not None:
        args.frac_random_rotamers_per_cst = [args.frac_random_rotamers]+[args.frac_random_rotamers for r in restypes]
    
        # In case best rotamer is requested for a given CST id then set randomness to 1.0
        for i, frac in enumerate(args.frac_random_rotamers_per_cst):
            if i in args.use_best_rotamer_cstids:
                args.frac_random_rotamers_per_cst[i] = 1.0
    
    
    #### PARSING SECONDARY STRUCTURE LENGTHS #####
    if args.N_len_per_cst is None:
        args.N_len_per_cst = [None]+[args.N_len for r in restypes]
    else:
        assert all([isinstance(x, int) for x in args.N_len_per_cst])
        args.N_len_per_cst = [None]+ args.N_len_per_cst
    
    if args.C_len_per_cst is None:
        args.C_len_per_cst = [None]+[args.C_len for r in restypes]
    else:
        assert all([isinstance(x, int) for x in args.C_len_per_cst])
        args.C_len_per_cst = [None]+ args.C_len_per_cst


    # Loading favored rotamers for each used residue type in each CST block
    # This allows different rotamer sets to be stored if same residue type should be
    # on different secondary structures in different CST blocks
    # TODO: could also consider enabling different probabilities for different CST's or AA's?  <-- partially done
    if args.secstruct_per_cst is None:
        args.secstruct_per_cst = [None]+[args.secstruct for r in restypes]
    else:
        assert all([x in "EH-" for x in args.secstruct_per_cst])
        args.secstruct_per_cst = [None]+ args.secstruct_per_cst
    return args


def parse_motif_input(motif_input, cst_atoms, restypes):
    motifs = {}
    for motif_txt in motif_input.split(","):
        motif_cst_no = int(motif_txt.split(":")[0])
        if motif_cst_no != 1:
            sys.exit("External motif not supported for not-first CST's right now.")
        motif_resno = int(motif_txt.split(":")[1])
        motif_fp = motif_txt.split(":")[2]
        motifs[motif_cst_no] = {"resno": motif_resno,
                                "pose": pyr.pose_from_file(motif_fp),
                                "fp": motif_fp,
                                "atoms": None}
        motif_resname = motifs[motif_cst_no]["pose"].residue(motif_resno).name3()
        assert motif_resname in restypes[motif_cst_no], f"{motif_resname} not found in {restypes}"

        # Finding the CST atoms for a given CST
        for sub_cst_block in cst_atoms[motif_cst_no]:
            for per_aa_cstset in sub_cst_block:
                if motif_resname in [aa.split("-")[0] for aa in per_aa_cstset.keys()]:
                    motif_resname_full = [aa for aa in per_aa_cstset.keys() if aa.split("-")[0]==motif_resname][0]
                    motifs[motif_cst_no]["atoms"] = per_aa_cstset[motif_resname_full]
        if motifs[motif_cst_no]["atoms"] is None:
            print(cst_atoms)
            sys.exit("Unable to find correct motif atoms based on the corresponding CST definition")
    return motifs


def parse_rotamer_subsampling(args, cst_atoms):
    chi_subsampling_levels = {}
    __xtrachi_cst_def = {}
    _extra_chi_definitions = {}
    if args.extra_chi is not None:
        # 1:2,2:2,3:1,4:1
        __xtrachi = args.extra_chi.split(",")
        _extra_chi_definitions = {int(x.split(":")[0]): int(x.split(":")[1]) for x in __xtrachi}
    
    elif args.extra_chi_per_cst is not None:
        # CSTNO-1:2,2:2,3:1,4:1 CSTNO2-1:1,2:1 
        __xtrachi_cst = {int(x.split("-")[0]): x.split("-")[1].split(",") for x in args.extra_chi_per_cst}
        __xtrachi_cst_def = {cstno: {int(x.split(":")[0]): int(x.split(":")[1]) for x in val} for cstno, val in __xtrachi_cst.items()}

    for cstno in cst_atoms:
        chi_subsampling_levels[cstno] = {}
        for n in range(4):
            if cstno in __xtrachi_cst_def.keys() and n+1 in __xtrachi_cst_def[cstno].keys():
                chi_subsampling_levels[cstno][n+1] = __xtrachi_cst_def[cstno][n+1]
            elif n+1 in _extra_chi_definitions.keys():
                chi_subsampling_levels[cstno][n+1] = _extra_chi_definitions[n+1]
            else:
                chi_subsampling_levels[cstno][n+1] = 0
            assert 0 <= chi_subsampling_levels[cstno][n+1] <= 7, f"Invalid sampling level for cst {cstno}, chi {n+1}: {chi_subsampling_levels[cstno][n+1]}"

    print("Using CHI sampling levels for CST's:")
    for cstno in chi_subsampling_levels:
        print(f"    CST {cstno} :: {chi_subsampling_levels[cstno]}")

    return chi_subsampling_levels


def get_cst_atoms(cst_io):
    cst_atoms = {}
    for n in range(1, cst_io.mcfi_lists_size()+1):
        cst_atoms[n] = []
        for m in range(1, cst_io.mcfi_list(n).num_mcfis()+1):
            cst_atoms[n].append([])
            _mcfi = cst_io.mcfi_list(n).mcfi(m)

            # Figuring out if there is a particular downstream or upstream secondary match happening
            downstream_match = False
            upstream_match = False
            downstream_res_cst = 1
            if _mcfi.algorithm_inputs().__contains__("match"):
                if any(["DOWNSTREAM" in ai for ai in _mcfi.algorithm_inputs()["match"]]):
                    downstream_match = True
                    downstream_res_cst = 1  # I think this is always 1, right?
                elif any(["UPSTREAM_CST" in ai for ai in _mcfi.algorithm_inputs()["match"]]):
                    upstream_match = True
                    for ai in _mcfi.algorithm_inputs()["match"]:
                        if "SECONDARY_MATCH:" in ai and "UPSTREAM_CST" in ai:
                            downstream_res_cst = int(ai.split()[2])
                            break


            rt_combs = itertools.product(_mcfi.allowed_restypes(_mcfi.downstream_res()), _mcfi.allowed_restypes(_mcfi.upstream_res()))
            for (ds_res, us_res) in rt_combs:
                ais_ds = [ds_res.atom_name(_mcfi.template_atom_inds(_mcfi.downstream_res(), ai, ds_res)[1]) for ai in range(1, 4)]
                ais_us = [us_res.atom_name(_mcfi.template_atom_inds(_mcfi.upstream_res(), ai, us_res)[1]) for ai in range(1, 4)]
                
                # Need to append CST numbers to residue names
                cst_atoms[n][-1].append({f"{ds_res.name()}-{downstream_res_cst}": tuple(ais_ds),
                                         f"{us_res.name()}-{n}": tuple(ais_us)})

    return cst_atoms



"""
ROTAMER-RELATED FUNCTIONS
"""
def preselect_inverse_rotamers(rotset, restype_good_rotamers, keep_his_tautomer_per_cst, tip_atom=False):
    if tip_atom is False:
        print("Preselecting inverse rotamers based on Dunbrack probability")
        good_rotamers = [[] for x in rotset]
        for i, invrots in enumerate(rotset):
            if len(invrots) == 0:
                continue
            for res in invrots:
                if isinstance(res, pyrosetta.rosetta.core.pose.Pose):  # motif pose
                    good_rotamers[i].append(res)
                    continue
                if res.is_ligand():
                    # if len(good_rotamers[i]) > 0 and args.single_ligand_rotamer is True:
                    #     break
                    good_rotamers[i].append(res)
                    continue
                if res.name3() == "HIS" and keep_his_tautomer_per_cst is not None:
                    if res.name() != keep_his_tautomer_per_cst[i]:
                        continue
                # Need to exclude proton CHIs
                _chis = [res.chi(n+1) for n in range(res.nchi()) if "H" not in [res.atom_type(an).element() for an in res.chi_atoms(n+1)]]
                if res.name3() in ["ALA", "GLY"]:
                    good_rotamers[i].append(res)
                else:
                    rotlib_matches = dunbrack_rotlib.find_bb_from_inverse_loc(restype_good_rotamers[i][res.name3()], _chis)
                    if len(rotlib_matches) > 0:
                        good_rotamers[i].append(res)
            if len(good_rotamers[i]) == 0 and len(rotset[i]) != 0:
                print(f"Failed to find compatible rotamers for constraint {i}: {res.name()}")
                return None
    else:
        print("Preselecting inverse rotamers only based whether the tip atoms are different")
        good_rotamers = []
        for i, invrots in enumerate(rotset):
            if isinstance(invrots[0], pyrosetta.rosetta.core.pose.Pose):
                good_rotamers.append(invrots)
                continue
            elif invrots[0].is_ligand():
                good_rotamers.append(invrots)
                continue
            good_rotamers.append([])
            for invrot in invrots:
                if len(good_rotamers[i]) == 0:
                    good_rotamers[i].append(invrot)
                    continue
                is_unique = []
                for rot in good_rotamers[i]:
                    if rot.name() != invrot.name():
                        continue
                    if (rot.xyz("CA")-invrot.xyz("CA")).norm() < 0.2:
                        is_unique.append(False)
                        continue
                    if (rot.xyz("CB")-invrot.xyz("CB")).norm() < 0.2:
                        is_unique.append(False)
                        continue
                    is_unique.append(True)
                if all(is_unique):
                    good_rotamers[i].append(invrot)
    return good_rotamers


def find_unique_rotamers_for_motif(rotset, motifs):
    """
    Identifies different rotamers from the inverse rotamer set that can be used for aligning the motif to.
    Difference is calculated based on the geometric distance between the motif atoms of inverse rotamers.
    """
    print("Preselecting inverse rotamers for motif alignment, based on unique CST subsampling")
    unique_rotset = []
    
    for i, invrots in enumerate(rotset):
        if len(invrots) == 0:
            unique_rotset.append([])
            continue
        unique_rotamers = []
        for j, res in enumerate(invrots):
            if len(unique_rotamers) == 0:
                unique_rotamers.append(res)
                continue
            dms = []
            for ures in unique_rotamers:
                dms.append([(res.xyz(a)-ures.xyz(a)).norm() for a in motifs[i]["atoms"]])
            if all([sum(x) > 0.1 for x in dms]):
                unique_rotamers.append(res)

        print(f"    CST {i}, {len(unique_rotamers)}/{len(invrots)} after unique selection")
        unique_rotset.append(unique_rotamers)
    return unique_rotset


def pick_random_rotamers(invrots, N_max=None, frac=None):
    if N_max is not None:
        if len(invrots) < N_max:
            return [r for r in invrots]
        elif isinstance(invrots[0], pyrosetta.rosetta.core.pose.Pose):
            return [r for r in invrots]
        else:
            return random.sample([r for r in invrots], N_max)
    if frac is not None:
        if len(invrots) <= 1:
            return [r for r in invrots]
        elif isinstance(invrots[0], pyrosetta.rosetta.core.pose.Pose):
            return [r for r in invrots]
        else:
            return random.sample([r for r in invrots], int(round(frac*len(invrots), 0)))


def pick_random_rotamers_set(rotset, max_random_rotamers_per_cst=None, frac_random_rotamers_per_cst=None):
    """
    Selects a subset of inverse rotamers for each set of inverse rotamers
    Arguments:
        rotset (list)
        max_random_rotamers_per_cst (list, int)
        frac_random_rotamers_per_cst (list, float)
    """
    if max_random_rotamers_per_cst is None and frac_random_rotamers_per_cst is None:
        sys.exit("Bad setup")
    elif max_random_rotamers_per_cst is not None and frac_random_rotamers_per_cst is not None:
        sys.exit("Bad setup")

    if max_random_rotamers_per_cst is None:
        max_random_rotamers_per_cst = [None for x in frac_random_rotamers_per_cst]
    elif frac_random_rotamers_per_cst is None:
        frac_random_rotamers_per_cst = [None for x in max_random_rotamers_per_cst]

    assert len(rotset) == len(frac_random_rotamers_per_cst)
    assert len(rotset) == len(max_random_rotamers_per_cst)

    rotsett = []

    for n, invrots in enumerate(rotset):
        rotsett.append(pick_random_rotamers(invrots, N_max=max_random_rotamers_per_cst[n], frac=frac_random_rotamers_per_cst[n]))
    return rotsett


def subsample_rotamers(rotamers, subsample_levels, per_cst_rotlib, cst_atoms):
    expanded_rotset = []
    for cst_block, invrots in enumerate(rotamers):
        expanded_rotset.append([])
        if cst_block == 0:  # Ligand
            expanded_rotset[0] = [r for r in invrots]
            continue
        for n, invrot in enumerate(invrots):
            if isinstance(invrot, pyrosetta.rosetta.core.pose.Pose):  # motif pose
                expanded_rotset[cst_block].append(invrot)
                continue
            _asd = dunbrack_rotlib.find_bb_from_inverse_loc(per_cst_rotlib[cst_block][invrot.name3()], list(invrot.chi()))
            if len(_asd) == 0:
                print(f"CST {cst_block}: rotamer {n} found no hits from Dunbrack library!?")
                expanded_rotset[cst_block].append(invrot)
                continue
            # Right not taking STDEV just as an average of all found rotamers in desired secondary structure bins
            stdevs = {chino+1: _asd[f"std{chino+1}"].mean() for chino in range(invrot.nchi())}

            # Expanding all chi's based on user request
            chi_samplings = {chino: calculate_samplings(invrot.chi(chino), stdevs[chino], subsample_levels[cst_block][chino]) for chino in stdevs}
            for chiset in itertools.product(*chi_samplings.values()):
                _rot = invrot.clone()
                for chino, _chi in enumerate(chiset):
                    _rot.set_chi(chino+1, _chi)
                
                # Need to realign coordinates
                # First let's find what are the CST atoms used
                align_atoms = [[restype_block[f"{invrot.name()}-{cst_block}"] for restype_block in var_cst if invrot.name() == list(restype_block.keys())[1].split("-")[0]] for var_cst in cst_atoms[cst_block]]
                align_atoms = list(set([item for sublist in align_atoms for item in sublist]))
                if len(align_atoms) != 1:
                    print(f"Bad choice for alignment atoms: {align_atoms}")
                __rot = align_pdbs.align_residue_to_residue(invrot, _rot, {"atoms1": align_atoms[0],
                                                                           "atoms2": align_atoms[0]})
                expanded_rotset[cst_block].append(__rot)
        print(f"Expanded CST-{cst_block} rotamers from {len(invrots)} to {len(expanded_rotset[cst_block])}")
    return expanded_rotset



def prune_ligand_rotamers(rotset, rmsd_cutoff=None, nproc=None):
    print("Pruning ligand rotamers based on intramolecular clashes")
    # Clashcheck
    def process():
        while True:
            i = the_queue.get(block=True)
            if i is None:
                return
            res = rotset[i]
            nonbonded_distmat = []
            for p in itertools.combinations(range(1, res.natoms()+1), 2):
                if any([res.is_virtual(n) for n in p]):
                    continue
                # Skipping over bonded atoms
                if p[0] in res.bonded_neighbor(p[1]) or p[1] in res.bonded_neighbor(p[0]):
                    continue
                nonbonded_distmat.append((res.xyz(p[0]) - res.xyz(p[1])).norm())

                if all([res.atom_type(n).is_heavyatom() for n in p]):
                    cutoff = 2.1
                else:
                    cutoff = 1.7

                if nonbonded_distmat[-1] < cutoff:
                    # if args.debug: print(f"Ligand rotamer pruning: {i}: {p}, {res.atom_name(p[0])}-{res.atom_name(p[1])}, {nonbonded_distmat[-1]}")
                    good_rotamers[i] = False
                    # print(f"Clashing ligand rotamer: {i}")
                    break

    print(f"{len(rotset)} conformers to process")
    the_queue = multiprocessing.Queue()  # Queue stores the iterables

    start = time.time()
    manager = multiprocessing.Manager() 
    good_rotamers = manager.dict()  # Need a special dictionary to store outputs from multiple processes

    for i, res in enumerate(rotset):
        the_queue.put(i)
        good_rotamers[i] = True

    pool = multiprocessing.Pool(processes=nproc,
                                initializer=process)

    # None to end each process
    for _i in range(nproc):
        the_queue.put(None)

    # Closing the queue and the pool
    the_queue.close()
    the_queue.join_thread()
    pool.close()
    pool.join()

    end = time.time()
    print(f"Found {len([i for i in good_rotamers.keys() if good_rotamers[i] is True])} good ligand rotamers.\n"
          f"Processing all the rotamers for clashes took {(end - start):.2f} seconds")
    
    
    ## RMSD
    if rmsd_cutoff in [None, 0.0]:
        return [rotset[i] for i in good_rotamers.keys() if good_rotamers[i] is True]

    unique_rotamers = {}
    DMs = {}
    for i in good_rotamers.keys():
        if good_rotamers[i] is False:
            continue
        res = rotset[i]

        xyz = np.array([res.xyz(n+1) for n in range(res.natoms()) if res.atom_type(n+1).element() != "H"])
        DMs[i] = scipy.spatial.distance.pdist(xyz, 'euclidean')

        if len(unique_rotamers) == 0:
            unique_rotamers[i] = res
            continue
        rmsds = []
        for j, res_u in unique_rotamers.items():
            rmsds.append(utils.rmsd(DMs[i], DMs[j]))
            
            if rmsds[-1] < rmsd_cutoff:
                break

        if min(rmsds) < rmsd_cutoff:
            continue
        else:
            unique_rotamers[i] = rotset[i]

    print(f"Found {len(unique_rotamers)}/{len(good_rotamers)} unique ligand rotamers based on RMSD cutoff {rmsd_cutoff}.")
    return [rot for i, rot in unique_rotamers.items()]


def prune_residue_rotamers(rotset):
    """
    Pruning based on proton chi similarity
    """
    unique_rotamers = {}
    for i, res in enumerate(rotset):
        if res.name3() not in utils.N_chis:
            n_chis = len([n for n in range(1, res.nchi()+1) if not any([res.atom_type(x).element() == "H" for x in res.chi_atoms(n)])])
        else:
            n_chis = utils.N_chis[res.name3()]
        if res.nchi() == n_chis:
            unique_rotamers[i] = res
            continue
        if i == 0:
            unique_rotamers[i] = res
            continue

        ads = []  # largest atom-atom distance between heavyatoms of RES and all parsed residues
        for j, res_u in unique_rotamers.items():
            if res.name3() != res_u.name3():
                continue
            ads.append(max([(res.xyz(n+1) - res_u.xyz(n+1)).norm() for n in range(res.natoms()) if res.atom_type(n+1).element() != "H"]))
            if ads[-1] < 0.02:
                break

        if len(ads) == 0 or min(ads) >= 0.02:
            unique_rotamers[i] = res
        else:
            continue

    return [val for k, val in unique_rotamers.items()]


def calculate_samplings(chi_value, std, sampling_level):
    """
    0 Default  original dihedral only; same as using no flag at all
    1          +/- one standard deviation (sd); 3 samples
    2          +/- 0.5 sd; 3 samples
    3          +/- 1 & 2 sd; 5 samples
    4          +/- 0.5 & 1 sd; 5 samples
    5          +/- 0.5, 1, 1.5 & 2 sd; 9 samples
    6          +/- 0.33, 0.67, 1 sd; 7 samples
    7          +/- 0.25, 0.5, 0.75, 1, 1.25 & 1.5 sd; 13 samples.
    """
    if sampling_level == 0:
        samples = [chi_value]
    elif sampling_level == 1:
        samples = [chi_value-std, chi_value, chi_value+std]
    elif sampling_level == 2:
        samples = [chi_value-0.5*std, chi_value, chi_value+0.5*std]
    elif sampling_level == 3:
        samples = [chi_value-2*std, chi_value-std, chi_value, chi_value+std, chi_value+2*std]
    elif sampling_level == 4:
        samples = [chi_value-std, chi_value-0.5*std, chi_value, chi_value+0.5*std, chi_value+std]
    elif sampling_level == 5:
        samples = [chi_value-2*std, chi_value-1.5*std, chi_value-std, chi_value-0.5*std, 
                   chi_value,
                   chi_value+0.5*std, chi_value+std, chi_value+1.5*std, chi_value+2*std]
    elif sampling_level == 6:
        samples = [chi_value*std, chi_value-0.667*std, chi_value-0.333*std, 
                   chi_value,
                   chi_value+0.333*std, chi_value+0.667*std, chi_value*std]
    elif sampling_level == 7:
        samples = [chi_value-1.5*std, chi_value-1.25*std, chi_value-std, chi_value-0.75*std, chi_value-0.5*std, chi_value-0.25*std,
                   chi_value,
                   chi_value+0.25*std, chi_value+0.5*std, chi_value+0.75*std, chi_value+std, chi_value+1.25*std, chi_value+1.5*std]
    else:
        sys.exit(f"Invalid sampling level: {sampling_level}")
    return samples


"""
Functions used during inverse rotamer assembly generation
"""

def check_clash(pose, catres_resnos, cutoff=1.7, ignore_respairs=None, tip_atom=False, debug=False):
    """
    Checks for clashes between residue atoms
    Only consideres residues that have nbr_atom within 10 angstrom of eachother.
    Default clash cutoff is 1.7 angstrom.
    Clashes are not detected for N-H and O-H contacts.
    """

    # combs = itertools.product(*[x for x in [pose.residues, pose.residues]])
    combs = itertools.combinations(range(1, pose.size()+1), 2)
    for c in combs:
        res1 = pose.residue(c[0])
        res2 = pose.residue(c[1])
        # Going through a bunch of conditions that would allow us to skip
        # checking clashes in a given pair of residues
        
        _ignore_atoms = {res1.seqpos():[],res2.seqpos():[]}
        if tip_atom is True:
            # Ignoring any of the backbone-ish atoms
            for r in [res1, res2]:
                if r.is_ligand():
                    continue
                if r.seqpos() in catres_resnos:
                    if r.name3() in ["GLY", "PRO", "ALA"]:
                        continue
                    for a in ["CA", "CB", "C", "N", "O"]:
                        _ignore_atoms[r.seqpos()].append(r.atom_index(a))
                        if r.attached_H_begin(r.atom_index(a)) == 0:
                            continue
                        for _n in range(r.attached_H_begin(r.atom_index(a)), r.attached_H_end(r.atom_index(a))+1):
                            _ignore_atoms[r.seqpos()].append(_n)

        if ignore_respairs is not None:
            if any([res1.seqpos() in p and res2.seqpos() in p for p in ignore_respairs]):
                continue

        if res1.chain() == res2.chain():
            continue
        if res1.seqpos() == res2.seqpos():
            continue
        if res1.is_bonded(res2):
            continue
        if (res1.nbr_atom_xyz() - res2.nbr_atom_xyz()).norm() > 10.0:
            continue
        if res1.is_virtual_residue() or res2.is_virtual_residue():
            continue
        for atm1 in range(1, res1.natoms()+1):
            if res1.is_virtual(atm1):
                continue
            if atm1 in  _ignore_atoms[res1.seqpos()]:
                continue
            for atm2 in range(1, res2.natoms()+1):
                if res2.is_virtual(atm2):
                    continue
                if atm2 in  _ignore_atoms[res2.seqpos()]:
                    continue

                if all([res1.atom_type(atm1).is_heavyatom(), res2.atom_type(atm2).is_heavyatom()]):
                    cutoff = 1.8
                else:
                    cutoff = 1.6
                _dist = (res1.xyz(atm1) - res2.xyz(atm2)).norm()
                if _dist < cutoff:
                    if res1.atom_type(atm1).element() in "NO" and res2.atom_type(atm2).element() == "H":  # H-bonds are not clashes
                        continue
                    # elif res1.atom_type(atm1).element() == "H" and res2.atom_type(atm2).element() in "NO":
                    #     continue
                    else:
                        if debug: print(f"Clashing atoms: {res1.name()}-{res1.seqpos()}-{res1.atom_name(atm1)} -- {res2.name()}-{res2.seqpos()}-{res2.atom_name(atm2)}: {_dist}")
                        return True
    return False


def adjust_bb(pose, resno, phi, psi):
    pose.set_phi(resno, phi)
    pose.set_psi(resno, psi)
    pose.set_omega(resno, 180.0)


def extend_SS(pose, ref_seqpos, secstruct, AAA, nres_Nterm=4, nres_Cterm=5):
    """
    Extends the stubs around a given residue in a pose by a number of residues on N-term and C-term side.
    The secondary structure is set to either idealized Helix or Strand

    Parameters
    ----------
    pose : pyrosetta.rosetta.core.pose.Pose
        DESCRIPTION.
    ref_seqpos : int
        DESCRIPTION.
    secstruct : str
        "E" or "H".
    AAA : pyrosetta.rosetta.core.pose.Pose
        pose object with 3 alanines.
    nres_Nterm : int, optional
        How many residues are added to N terminus. The default is 4.
    nres_Cterm : int, optional
        How many residues are added to C terminus. The default is 5.

    Returns
    -------
    pose2 : TYPE
        DESCRIPTION.

    """
    # assert nres_Nterm >= 2, "Too short N-term extension"
    # assert nres_Cterm >= 2, "Too short C-term extension"
    pose2 = pose.clone()
    for n in range(nres_Cterm):
        pose2.append_polymer_residue_after_seqpos(AAA.residue(2), ref_seqpos+n, True)
        adjust_bb(pose2, ref_seqpos+n, phi=utils.idealized_SS_phi_psi[secstruct]["phi"][0], psi=utils.idealized_SS_phi_psi[secstruct]["psi"][0])

    if nres_Cterm > 0:
        adjust_bb(pose2, pose2.size(), phi=utils.idealized_SS_phi_psi[secstruct]["phi"][0], psi=utils.idealized_SS_phi_psi[secstruct]["psi"][0])
    else:
        # If no C-term stub included then adding temporarily one
        pose2.append_polymer_residue_after_seqpos(AAA.residue(2), ref_seqpos, True)

    for n in range(nres_Nterm):
        pose2.prepend_polymer_residue_before_seqpos(AAA.residue(2), ref_seqpos, True)
        if n == 0:
            # Building foldtree to have a center point at the reference residue
            ft = pyrosetta.rosetta.core.kinematics.FoldTree()
            ft.add_edge(ref_seqpos+2, pose2.chain_begin(pose2.chain(ref_seqpos)), -1)
            ft.add_edge(ref_seqpos+2, pose2.chain_end(pose2.chain(ref_seqpos)), -1)
            for j in range(1, pose2.num_chains()+1):
                if j == pose2.chain(ref_seqpos):
                    continue
                else:  # adding foldtree edges for other chains
                    ft.add_edge(pose2.fold_tree().get_residue_edge(pose2.chain_begin(j)))
            pose2.fold_tree().clear()
            pose2.fold_tree(ft)
        adjust_bb(pose2, ref_seqpos+1, phi=utils.idealized_SS_phi_psi[secstruct]["phi"][0], psi=utils.idealized_SS_phi_psi[secstruct]["psi"][0])

    adjust_bb(pose2, ref_seqpos, phi=utils.idealized_SS_phi_psi[secstruct]["phi"][0], psi=utils.idealized_SS_phi_psi[secstruct]["psi"][0])

    if nres_Cterm == 0:
        pose2.delete_residue_slow(pose2.size())

    return pose2


def create_remark_lines(pose, catalytic_residues, cst_io):
    ## Adding REMARK 666 lines to the PDB's
    ## This is actually quite arduous since we need to figure out which variable CST block a particular residue came from

    pdb_info = pyrosetta.rosetta.core.pose.PDBInfo(pose)  # can this be added to pose somehow?

    ligands = [r for r in pose.residues if r.is_ligand()]

    calculators = {"dis": utils.get_dist, "ang": utils.get_angle, "tor": utils.get_dihedral}
    remarks = []
    for j, resno in catalytic_residues.items():
        if pose.residue(resno).is_ligand() and j == 0:
            continue
        rmrk = None

        for m in range(1, cst_io.mcfi_list(j).num_mcfis()+1):
            _mcfi = cst_io.mcfi_list(j).mcfi(m)

            downstream_res_cst = 0
            if _mcfi.algorithm_inputs().__contains__("match"):
                if any(["DOWNSTREAM" in ai for ai in _mcfi.algorithm_inputs()["match"]]):
                    downstream_res_cst = 0  # I think this is always 1, right?
                elif any(["UPSTREAM_CST" in ai for ai in _mcfi.algorithm_inputs()["match"]]):
                    for ai in _mcfi.algorithm_inputs()["match"]:
                        if "SECONDARY_MATCH:" in ai and "UPSTREAM_CST" in ai:
                            downstream_res_cst = int(ai.split()[2])
                            break
            # Residues in the final pose
            DS_RES = pose.residue(catalytic_residues[downstream_res_cst])
            US_RES = pose.residue(resno)
            
            good_cst_found = False
            rt_combs = itertools.product(_mcfi.allowed_restypes(_mcfi.downstream_res()), _mcfi.allowed_restypes(_mcfi.upstream_res()))
            for (ds_res, us_res) in rt_combs:
                if US_RES.name().split(":")[0] != us_res.name():  # skipping the wrong residue types
                    continue
                ais_ds = [ds_res.atom_name(_mcfi.template_atom_inds(_mcfi.downstream_res(), ai, ds_res)[1]) for ai in range(1, 4)]
                ais_us = [us_res.atom_name(_mcfi.template_atom_inds(_mcfi.upstream_res(), ai, us_res)[1]) for ai in range(1, 4)]

                cst_atomsets = {'dis_U1D1': [DS_RES.xyz(ais_ds[0]), US_RES.xyz(ais_us[0])],
                                'ang_U1D2': [DS_RES.xyz(ais_ds[1]), DS_RES.xyz(ais_ds[0]), US_RES.xyz(ais_us[0])],
                                'ang_U2D1': [DS_RES.xyz(ais_ds[0]), US_RES.xyz(ais_us[0]), US_RES.xyz(ais_us[1])],
                                'tor_U1D3': [DS_RES.xyz(ais_ds[2]), DS_RES.xyz(ais_ds[1]), DS_RES.xyz(ais_ds[0]), US_RES.xyz(ais_us[0])],
                                'tor_U2D2': [DS_RES.xyz(ais_ds[1]), DS_RES.xyz(ais_ds[0]), US_RES.xyz(ais_us[0]), US_RES.xyz(ais_us[1])],
                                'tor_U3D1': [DS_RES.xyz(ais_ds[0]), US_RES.xyz(ais_us[0]), US_RES.xyz(ais_us[1]), US_RES.xyz(ais_us[2])]}
                cst_atomsets = {k: np.array(v) for k,v in cst_atomsets.items()}
                
                # Measuring whether a particular respair geometrically matches the CST
                good_cst_found = False
                for cs in _mcfi.constraints():
                    passed_cst = []
                    for cst_par in cst_atomsets.keys():
                        cst_samples = getattr(cs, cst_par).create_sample_vector()
                        val = calculators[cst_par[:3]](*cst_atomsets[cst_par])
                        if val < 0.0:
                            val = 360.0 + val
                        # is any of the sampled values very close to the measured value?
                        if "dis" in cst_par:
                            passed_cst.append( any([abs(val-x) < 0.1 for x in cst_samples]) )
                        else:
                            passed_cst.append( any([abs(val-x) < 1.0 for x in cst_samples]) )
                    if all(passed_cst):
                        good_cst_found = True
                        break
                if good_cst_found:
                    break
            if good_cst_found:
                # if there's only one ligand then it will be stored as chain X residue 0
                if len(ligands) == 1 and DS_RES.name3() == ligands[0].name3():
                    rmrk = f"REMARK 666 MATCH TEMPLATE X {DS_RES.name3()}"\
                            f"    0 MATCH MOTIF {pdb_info.chain(resno)} "\
                            f"{US_RES.name3()} {resno:>4}  {j}  {m}               "
                else:
                    rmrk = f"REMARK 666 MATCH TEMPLATE {pdb_info.chain(DS_RES.seqpos())} {DS_RES.name3()}"\
                            f" {DS_RES.seqpos():>4} MATCH MOTIF {pdb_info.chain(resno)} "\
                            f"{US_RES.name3()} {resno:>4}  {j}  {m}               "
                remarks.append(rmrk)
                break
    return remarks

