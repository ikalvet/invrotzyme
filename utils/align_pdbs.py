import os,sys 
sys.path.append(os.path.dirname(os.path.realpath(__file__)))
import kabsch_align 
import util 
import numpy as np
import pyrosetta as pyr
import pyrosetta.rosetta


def find_atom_idx(atom, mapping):
    for i,A in enumerate(mapping):
        try:
            if A.strip() == atom:
                return i
        except AttributeError:
            print('This is atom ',A)

    raise KeyError(f'Could not find atom {atom} in mapping {mapping}')


def align_pose_to_residue(ref_residue, mobile_pose, ref_atoms):
    xyz1, parsed1 = get_xyz_stack_residue(ref_residue, ref_atoms["atoms1"])
    xyz2, parsed2 = get_xyz_stack_pose(mobile_pose, ref_atoms["atoms2"])

    # run Kabsch to get rotation matrix for atoms and rmsd
    # aligns xyz2 onto xyz1
    rmsd, _, R = kabsch_align.np_kabsch(xyz1, xyz2)
    print('RMSD between atoms: ',rmsd)

    # (1) now translate both proteins such that centroid(xyz1/xyz2) is at origin
    # (2) rorate xyz2 onto xyz1 with R
    # (3) write pdbs into outdir

    def centroid(X):
        # return the mean X,Y,Z down the atoms
        return np.mean(X, axis=0, keepdims=True)

    # centroid of just the points being aligned
    centroid1 = centroid(xyz1)
    centroid2 = centroid(xyz2)

    # (1)
    #xyz_protein1 = np.copy(parsed1['xyz']) - centroid1
    xyz_protein2 = np.copy(parsed2) - centroid2

    # (2)
    xyz_protein2 = xyz_protein2 @ R

    # Translate protein 2 to where it aligns with original protein 1
    xyz_protein2 += centroid1
    
    out_pose = mobile_pose.clone()
    for resno, res_coords in enumerate(xyz_protein2):
        for i, ac in enumerate(res_coords):
            if np.isnan(ac[0]):
                break
            out_pose.residue(resno+1).set_xyz(i+1, pyrosetta.rosetta.numeric.xyzVector_double_t(*ac))
            continue
    return out_pose


def align_residue_to_residue(ref_residue, mobile_residue, ref_atoms):
    xyz1, parsed1 = get_xyz_stack_residue(ref_residue, ref_atoms["atoms1"])
    xyz2, parsed2 = get_xyz_stack_residue(mobile_residue, ref_atoms["atoms2"])

    # run Kabsch to get rotation matrix for atoms and rmsd
    # aligns xyz2 onto xyz1
    rmsd, _, R = kabsch_align.np_kabsch(xyz1, xyz2)
    if rmsd > 0.1:
        print('RMSD between atoms: ',rmsd)

    # (1) now translate both proteins such that centroid(xyz1/xyz2) is at origin
    # (2) rorate xyz2 onto xyz1 with R
    # (3) write pdbs into outdir

    def centroid(X):
        # return the mean X,Y,Z down the atoms
        return np.mean(X, axis=0, keepdims=True)

    # centroid of just the points being aligned
    centroid1 = centroid(xyz1)
    centroid2 = centroid(xyz2)

    # (1)
    #xyz_protein1 = np.copy(parsed1['xyz']) - centroid1
    xyz_protein2 = np.copy(parsed2) - centroid2

    # (2)
    xyz_protein2 = xyz_protein2 @ R

    # Translate protein 2 to where it aligns with original protein 1
    xyz_protein2 += centroid1
    
    out_residue = mobile_residue.clone()

    for i, ac in enumerate(xyz_protein2[0]):
        if np.isnan(ac[0]):
            break
        out_residue.set_xyz(i+1, pyrosetta.rosetta.numeric.xyzVector_double_t(*ac))
        continue
    return out_residue


def get_xyz_stack_residue(residue, atoms_list):
    """
    Extracts the xyz crds corresponding to every atom in atoms_list 
    atoms_list format: [(resno, atomname), (resno, atomname), ...]
    """
    if residue.is_ligand() or residue.is_virtual_residue():
        return None, None

    xyz_all = parse_residue_coords(residue)
    seq = [util.alpha_1.index(residue.name1())]
    xyz_out = []

    # for each atom, get residue index and atom index 
    # store crds 
    for atom in atoms_list:
        # get index of residue and its Heavy atom mapping
        AA_int = seq[0]

        if residue.is_lower_terminus():
            AA_long_map = util.aa2longH_Nterm[AA_int]
        elif residue.is_upper_terminus():
            AA_long_map = util.aa2longH_Cterm[AA_int]
        else:
            AA_long_map = util.aa2longH[AA_int]

        # get index of atom in residue 
        atom_idx0 = find_atom_idx(atom.strip(), AA_long_map)

        # crds of this atom 
        xyz_atom = xyz_all[0, atom_idx0, :]

        xyz_out.append(xyz_atom)

    return np.array(xyz_out), xyz_all


def get_xyz_stack_pose(pose, atoms_list):
    """
    Extracts the xyz crds corresponding to every atom in atoms_list 
    atoms_list format: [(resno, atomname), (resno, atomname), ...]
    """

    xyz_all = parse_pose_coords(pose)
    seq = [util.alpha_1.index(r.name1()) for r in pose.residues if not r.is_ligand() and not r.is_virtual_residue()]
    xyz_out = []

    # for each atom, get residue index and atom index 
    # store crds 
    for (resn, atom) in atoms_list:
        # get index of residue and its Heavy atom mapping
        AA_int = seq[resn-1]
        if pose.residue(resn).is_lower_terminus():
            AA_long_map = util.aa2longH_Nterm[AA_int]
        elif pose.residue(resn).is_upper_terminus():
            AA_long_map = util.aa2longH_Cterm[AA_int]
        else:
            AA_long_map = util.aa2longH[AA_int]

        # get index of atom in residue 
        atom_idx0 = find_atom_idx(atom.strip(), AA_long_map)

        # crds of this atom 
        xyz_atom = xyz_all[resn-1, atom_idx0, :]

        xyz_out.append(xyz_atom)

    return np.array(xyz_out), xyz_all


def parse_pose_coords(pose):
    res = [r.seqpos() for r in pose.residues if not r.is_ligand() and not r.is_virtual_residue()]
    xyz = np.full((len(res), 26, 3), np.nan, dtype=np.float32)
    for r in pose.residues:
        if r.is_ligand() or r.is_virtual_residue():
            continue
        # rc = np.ndarray((res.natoms(), 3), dtype=np.float32)
        for n in range(r.natoms()):
            try:
                xyz[r.seqpos()-1][n] = r.xyz(n+1)
            except IndexError:
                print(r.name())
                print(r.seqpos())
                print(r.natoms())
                sys.exit(1)
    return xyz


def parse_residue_coords(residue):
    xyz = np.full((1, 26, 3), np.nan, dtype=np.float32)
    if residue.is_ligand() or residue.is_virtual_residue():
        return None
    # rc = np.ndarray((res.natoms(), 3), dtype=np.float32)
    for n in range(residue.natoms()):
        xyz[0][n] = residue.xyz(n+1)
    return xyz


