#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 25 23:12:52 2024

@author: indrek
"""
import numpy as np


# number of chis, excluding proton-chis
N_chis = {'ALA': 0, 'ARG': 4, 'TRP': 2, 'GLY': 0, 'ASP': 2, 'HIS': 2, 'GLU': 3,
          'GLN': 3, 'ASN': 2, 'LEU': 2, 'ILE': 2, 'THR': 1, 'VAL': 1, 'SER': 1,
          'MET': 3, 'CYS': 1, 'PRO': 3, 'LYS': 4, 'PHE': 2, 'TYR': 2, "CYX": 1}


# PHI and PSI values for ideal backbone, and tolerances for randomization
idealized_SS_phi_psi = {"H": {"phi": (-57.0, 10.0), "psi": (-47.0, 10.0)},
                        "E": {"phi": (-140.0, 20.0), "psi": (130.0, 20.0)},
                        "-": {"phi": (-140.0, 20.0), "psi": (130.0, 20.0)}}


def get_dist(a, b):
    return np.linalg.norm(a-b)


def get_angle(a1, a2, a3):
    a1 = np.array(a1)
    a2 = np.array(a2)
    a3 = np.array(a3)

    ba = a1 - a2
    bc = a3 - a2

    cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
    angle = np.arccos(cosine_angle)

    return round(np.degrees(angle), 1)



def get_dihedral(a1, a2, a3, a4):
    """
    a1, a2, a3, a4 (np.array)
    Each array has to contain 3 floats corresponding to X, Y and Z of an atom.
    Solution by 'Praxeolitic' from Stackoverflow:
    https://stackoverflow.com/questions/20305272/dihedral-torsion-angle-from-four-points-in-cartesian-coordinates-in-python#
    1 sqrt, 1 cross product
    Calculates the dihedral/torsion between atoms a1, a2, a3 and a4
    Output is in degrees
    """

    b0 = a1 - a2
    b1 = a3 - a2
    b2 = a4 - a3

    # normalize b1 so that it does not influence magnitude of vector
    # rejections that come next
    b1 /= np.linalg.norm(b1)

    # vector rejections
    # v = projection of b0 onto plane perpendicular to b1
    #   = b0 minus component that aligns with b1
    # w = projection of b2 onto plane perpendicular to b1
    #   = b2 minus component that aligns with b1
    v = b0 - np.dot(b0, b1)*b1
    w = b2 - np.dot(b2, b1)*b1

    # angle between v and w in a plane is the torsion angle
    # v and w may not be normalized but that's fine since tan is y/x
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return np.degrees(np.arctan2(y, x))


def rmsd(geom, target):
    return np.sqrt(((geom - target) ** 2).mean())


