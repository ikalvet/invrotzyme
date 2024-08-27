#!/usr/bin/env python3
import pandas as pd
import os

comparisons = {'<=': '__le__',
               '<': '__lt__',
               '>': '__gt__',
               '>=': '__ge__',
               '=': '__eq__'}

chi_psi_SS = {"H": {"phi": (-72.0, -50.0),
                    "psi": (-50.0, -30.0)},
              "E": {"phi": (-161.0, -89.0),
                    "psi": (109.0, 151.0)},
              "L": {"phi": (),
                    "psi": ()},
              "-": {"phi": (-180.0, 180.0),
                    "psi": (-180.0, 180.0)}}


def load_rotamer_df(dunbrack_database):
    header = ["restype", "phi", "psi", "N", "r1", "r2", "r3", "r4", "prob", "chi1", "chi2", "chi3", "chi4", "std1", "std2", "std3", "std4"]
    rotlib = pd.read_csv(dunbrack_database, sep="\s+", names=header)
    for n in range(1, 5):
        rotlib[f"chi{n}_min"] = rotlib[f"chi{n}"]-rotlib[f"std{n}"]
        rotlib[f"chi{n}_max"] = rotlib[f"chi{n}"]+rotlib[f"std{n}"]
    return rotlib


def filter_rotlib(scores, filters):
    filtered_scores = scores.copy()

    for s in filters.keys():
        _fltrs = []
        if isinstance(filters[s][0], list):
            _fltrs = filters[s]
        else:
            _fltrs.append(filters[s])
        for fltr in _fltrs:
            if fltr is not None and s in scores.keys():
                val = fltr[0]
                sign = comparisons[fltr[1]]
                filtered_scores =\
                  filtered_scores.loc[(filtered_scores[s].__getattribute__(sign)(val))]
    return filtered_scores


def find_good_rotamers(rotlib, restype, cumulative_prob=1.0, secstruct=None, phi=None, psi=None, keep_only_best=False):
    """
    Arugments:
        rotlib (pandas.DataFrame)
        restype (str) :: name3 of an amino acid in the rotamer library
        cumulative_prob (float) :: cumulative probability up to which rotamers are returned
        secstruct (str, ('H', 'E')) :: secondary structure type for which rotamers are searched.
        phi (tuple, (float, float)) :: min and max phi value for defining a subset of the library
        psi (tuple, (float, float)) :: min and max psi value for defining a subset of the library
        keep_only_best (bool) :: only the highest probability rotamer is returned for each phi/psi bin
    """
    assert isinstance(phi, (tuple, type(None)))
    assert isinstance(psi, (tuple, type(None)))
    assert secstruct in ("H", "E", "-", None), "Not implemented for other secondary structures yet"
    # assert restype not in ["ALA", "GLY"], "No rotamer library for ALA and GLY"
    assert not all([x is None for x in [secstruct, phi]]), "Must provide either secstruct letter OR phi and psi values"
    assert not all([x is None for x in [secstruct, psi]]), "Must provide either secstruct letter OR phi and psi values"

    if secstruct is not None:
        phi_limits = chi_psi_SS[secstruct]["phi"]
        psi_limits = chi_psi_SS[secstruct]["psi"]
    elif phi is not None and psi is not None:
        phi_limits = phi
        psi_limits = psi
    else:
        print("Both phi and psi need to be defined")
        return None

    filters = {'restype': [restype, '='],
               'phi': [[phi_limits[0], '>='], [phi_limits[1], '<=']],
               'psi': [[psi_limits[0], '>='], [psi_limits[1], '<=']]}

    SS_rotlib = filter_rotlib(rotlib, filters)
    phi_psi_bins = list(set([(row.phi, row.psi) for idx, row in SS_rotlib.iterrows()]))
    df = pd.DataFrame()
    for phi_psi_bin in phi_psi_bins:
        _df = SS_rotlib.loc[(SS_rotlib["phi"] == phi_psi_bin[0]) & (SS_rotlib["psi"] == phi_psi_bin[1])]
        if keep_only_best is True:
            _df2 = _df.iloc[0]
        else:
            if cumulative_prob == 1.0:
                _df2 = _df.copy()
            else:
                _df2 = _df.loc[_df.prob.cumsum() <= cumulative_prob]

                # Also adding the next most probable rotamer that would push the cumulative sum over the cutoff
                # This fixes the issue where no rotamers are returned when the cutoff is lower than the prob of the most likely rotamer
                if len(_df2) == 0:
                    idx_to_add = 0
                elif len(_df2) < len(_df):
                    idx_to_add = len(_df2)
                else:
                    idx_to_add = None
                if idx_to_add is not None:
                    _df2 = pd.concat([_df2, _df.iloc[idx_to_add]], ignore_index=True)
                    # _df2 = _df2.append(_df.iloc[idx_to_add], ignore_index=True)
        # df = df.append(_df2, ignore_index=True)
        df = pd.concat([df, _df2], ignore_index=True)
    return df


def find_bb_from_inverse(rotlib, chis):
    df = pd.DataFrame()
    for idx, row in rotlib.iterrows():
        _chi_matches = []
        for i, ch in enumerate(chis):
            _chi_matches.append(row[f"chi{i+1}"]-row[f"std{i+1}"] <= ch <= row[f"chi{i+1}"]+row[f"std{i+1}"])
        if all(_chi_matches):
            # df = df.append(row)
            df = pd.concat([df, row])
    return df


def find_bb_from_inverse_loc(rotlib, chis):
    """
    Finds
    Arguments:
        rotlib (pandas.DataFrame) :: rotamer library. Preferrably for a given amino acid.
        chis (list) :: list of chi values
    """
    assert isinstance(rotlib, pd.DataFrame)
    rl = rotlib.copy()
    for i, ch in enumerate(chis):
        rl = rl.loc[(rl[f"chi{i+1}_min"] <= ch) & (rl[f"chi{i+1}_max"] >= ch)]
    return rl



