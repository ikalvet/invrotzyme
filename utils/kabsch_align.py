import numpy as np
import copy
#Gyu Rie Lee
#Borrowed kabsch code and modified slightly for superimposition

#Use kabsch algorithm to align van der Mers with mainchain atoms (or given subset of coord)
#get transformation matrix from xyz1 and xyz2 (could be N-CA-C of residues)
#then use this to align residue+functional group
#xyz1/coord_for_align1 would be the reference
#IMPORTANT: xyz1_in is being copied inside as xyz1 because xyz1_in will be used repeatedly outside of this code


def np_kabsch(A,B):
    """
    Numpy version of kabsch algorithm. Superimposes B onto A

    Parameters:
        (A,B) np.array - shape (N,3) arrays of xyz crds of points


    Returns:
        rms - rmsd between A and B
        R - rotation matrix to superimpose B onto A
        rB - the rotated B coordinates
    """
    A = np.copy(A)
    B = np.copy(B)

    def centroid(X):
        # return the mean X,Y,Z down the atoms
        return np.mean(X, axis=0, keepdims=True)

    def rmsd(V,W, eps=1e-6):
        # First sum down atoms, then sum down xyz
        N = V.shape[-2]
        return np.sqrt(np.sum((V-W)*(V-W), axis=(-2,-1)) / N + eps)


    N, ndim = A.shape

    # move to centroid
    A = A - centroid(A)
    B = B - centroid(B)

    # computation of the covariance matrix
    C = np.matmul(A.T, B)

    # compute optimal rotation matrix using SVD
    U,S,Vt = np.linalg.svd(C)


    # ensure right handed coordinate system
    d = np.eye(3)
    d[-1,-1] = np.sign(np.linalg.det(Vt.T@U.T))

    # construct rotation matrix
    R = Vt.T@d@U.T

    # get rotated coords
    rB = B@R

    # calculate rmsd
    rms = rmsd(A,rB)

    return rms, rB, R


def kabsch_align_coords(xyz1, xyz2_in, mobile_coord):

#    xyz1 = copy.deepcopy(xyz1_in)
    xyz2 = copy.deepcopy(xyz2_in)
    # check dimensions
    #print(len(xyz1), len(xyz2))
    assert len(xyz1) == len(xyz2)
    L = len(xyz1)
    assert L > 2

    # move two both sets of points to their
    # centers of masses (COM)
    COM1 = np.sum(xyz1, axis=0) / float(L)
    COM2 = np.sum(xyz2, axis=0) / float(L)
    xyz1 -= COM1
    xyz2 -= COM2

    # Initial residual, see Kabsch.
    E0 = np.sum( np.sum(xyz1*xyz1,axis=0),axis=0) + np.sum( np.sum(xyz2*xyz2,axis=0),axis=0 )

    # SVD of the covariance matrix
    V, S, Wt = np.linalg.svd( np.dot( np.transpose(xyz2), xyz1))

    # check parity of the transformation
    reflect = float(str(float(np.linalg.det(V) * np.linalg.det(Wt))))
    if reflect == -1.0:
        S[-1] = -S[-1]
        V[:,-1] = -V[:,-1]

    RMSD = E0 - (2.0 * sum(S))
    RMSD = np.sqrt(abs(RMSD / L))

    # U is simply V*Wt
    U = np.dot(V, Wt)

    # translation vector
    t = COM1 - COM2

    superimposed_coord = np.dot((mobile_coord-COM2), U)
    superimposed_coord += COM1
#    rot_coord_2 = np.dot((coord_for_align2 - COM2), U)
#    rot_coord_1 = coord_for_align1 - COM1
    
#    rot_coord_2 = np.dot((coord_for_align2 - COM2), U) + COM1
    
#    return coord_for_align1, rot_coord_2
    return superimposed_coord
#    return RMSD, t, U

def kabsch_rmsd(xyz1_in,xyz2_in):

    xyz1 = copy.deepcopy(xyz1_in)
    xyz2 = copy.deepcopy(xyz2_in)
    # check dimensions
    assert len(xyz1) == len(xyz2)
    L = len(xyz1)
    assert L > 2

    # move two both sets of points to their
    # centers of masses (COM)
    COM1 = np.sum(xyz1, axis=0) / float(L)
    COM2 = np.sum(xyz2, axis=0) / float(L)
    xyz1 -= COM1
    xyz2 -= COM2

    # Initial residual, see Kabsch.
    E0 = np.sum( np.sum(xyz1*xyz1,axis=0),axis=0) + np.sum( np.sum(xyz2*xyz2,axis=0),axis=0 )

    # SVD of the covariance matrix
    V, S, Wt = np.linalg.svd( np.dot( np.transpose(xyz2), xyz1))

    # check parity of the transformation
    reflect = float(str(float(np.linalg.det(V) * np.linalg.det(Wt))))
    if reflect == -1.0:
        S[-1] = -S[-1]
        V[:,-1] = -V[:,-1]

    RMSD = E0 - (2.0 * sum(S))
    RMSD = np.sqrt(abs(RMSD / L))

    # U is simply V*Wt
    U = np.dot(V, Wt)

    # translation vector
    t = COM1 - COM2

    return RMSD
#    return RMSD, t, U

