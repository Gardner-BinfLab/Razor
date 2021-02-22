#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 17 10:36:23 2020
@author: bikash
"""

import re
import warnings
import numpy as np
import pandas as pd
from scipy.signal import savgol_filter

# Constants
S_MODELS = pd.read_pickle("libs/S.pkl.gz")
C_MODELS = pd.read_pickle("libs/C.pkl.gz")
FUNGI = pd.read_pickle("libs/Fungi_Classifier.pkl.gz")
TOXIN = pd.read_pickle("libs/Toxin_Classifier.pkl.gz")
weights_df = pd.read_pickle("libs/Cleavage_weights.pkl.gz")
WEIGHTS = [w.to_dict() for w in weights_df.Weight]
# Weights dataframe is of order weight['Position']['AA]
# where position is 0 based.

# Kyte & Doolittle index of hydrophobicity
# J. Mol. Biol. 157:105-132(1982).

# Normalized flexibility parameters (B-values),
# Vihinen M., Torkkila E., Riikonen P. Proteins. 19(2):141-9(1994).

# Solubility Weighted Index,
# Bhandari, B.K., Gardner, P.P. and Lim, C.S.,(2020),
# doi: 10.1093/bioinformatics/btaa578

hydrop_flex_swi = {
    "R": [-4.5, 1.008, 0.771],
    "K": [-3.9, 1.102, 0.927],
    "N": [-3.5, 1.048, 0.86],
    "D": [-3.5, 1.068, 0.908],
    "Q": [-3.5, 1.037, 0.789],
    "E": [-3.5, 1.094, 0.988],
    "H": [-3.2, 0.95, 0.895],
    "P": [-1.6, 1.049, 0.824],
    "Y": [-1.3, 0.929, 0.611],
    "W": [-0.9, 0.904, 0.637],
    "S": [-0.8, 1.046, 0.744],
    "T": [-0.7, 0.997, 0.81],
    "G": [-0.4, 1.031, 0.8],
    "A": [1.8, 0.984, 0.836],
    "M": [1.9, 0.952, 0.63],
    "C": [2.5, 0.906, 0.521],
    "F": [2.8, 0.915, 0.585],
    "L": [3.8, 0.935, 0.655],
    "V": [4.2, 0.931, 0.736],
    "I": [4.5, 0.927, 0.678],
}

def validate(seq, max_scan=45):
    """
    - Replaces 'U' with 'C'
    - Pads shorter sequence with 'S' so that the length
    is 30 residues.
    - Raises ValueError if 'X' is within the residues 
    defined by max_scan + 15.
    """
    seq = seq.upper()[: max_scan + 15].replace("U", "C")
    valid_aa = re.compile("^[RKNDQEHPYWSTGAMCFLVI]*$")
    match = re.match(valid_aa, seq)

    if match:
        length = len(seq)
        if length < 30:
            seq = seq + "S" * (30 - length)
        return seq
    else:
        raise ValueError(
            "Unknown residues in the input "
            "sequence.\n Only standard amino acid codes "
            "are allowed."
        )

def features(seq):
    """
    Features.
    Used to compute S score. So the sequence length should be 30.
    """

    seq = seq[:30]
    if len(seq) != 30:
        raise ValueError(
            "Input sequence must be 30 residues long!"
            "\nExpected length 30: Got {}".format(len(seq))
        )
    aa_list = 'RKNDCEVIYFWL' + 'QP'
    converted = np.array([hydrop_flex_swi[i] for i in seq])
    hydro = converted[:, 0]
    flex = converted[:, 1]
    swi = converted[:, 2]
    aa_counts = [seq.count(i) for i in aa_list]

    return np.concatenate(
        [savgol_filter(hydro, 15, 2), savgol_filter(swi, 15, 2), flex, aa_counts]
    )

def s_score(feat):
    """
    S score of sequence.
    Input is an array of features (104)
    """
    if len(feat) != 104:
        raise ValueError(
            "Input features length is incorrect!"
            "Expected length 104: Got {}".format(len(feat))
        )

    if feat.dtype != np.float64:
        raise TypeError("Non numeric characters not allowed!")

    classifiers = S_MODELS.Classifier
    S = np.array([clf.predict_proba([feat]) for clf in classifiers])[:, :, 1].flatten()
    return S

def validate_scan(seq, max_scan):

    if not type(max_scan) is int:
        raise TypeError("Only integers allowed for scan length.")
    if max_scan < 16:
        warnings.warn(
            "The minimum length to take for evaluating C score "
            "must be greater than 16 but received {max_scan}."
            " Correcting it to 45.".format(max_scan=max_scan)
        )
        max_scan = 45
    if max_scan > len(seq):
        warnings.warn(
            "The given maximum length to take for evaluating C score {max_scan} "
            "is greater than the input sequence length {len_seq}."
            " Correcting it to sequence length {len_seq}.".format(
                max_scan=max_scan, len_seq=len(seq)
            )
        )
        max_scan = len(seq)
    return max_scan

def c_score(seq, max_scan=45):
    """
    C score of sequence (Max probs in cleavage sites)
    Also returns the possible cleavage site and a probability of 
    cleavage sites along the sequence as scored by 5 models, 
    possible cleavage site (sites with max probs).
    """
    max_scan = validate_scan(seq, max_scan)

    if len(seq) <= len(seq[: max_scan + 15]):
        seq = seq + "S" * (15 - abs(len(seq) - max_scan))

    subseqs = [seq[i : i + 30] for i in range(max_scan - 15)]

    # Each subsequence is scored using each weight matrix.
    scored_subseqs = np.array(
        [
            [[weight[p][q] for p, q in enumerate(x)] for x in subseqs]
            for weight in WEIGHTS
        ]
    )

    # Each score is then classified by classifier corresponding to the weight.
    classifiers = C_MODELS.Classifier
    all_probs_ = np.array(
        [i[0].predict_proba(i[1]) for i in zip(classifiers, scored_subseqs)]
    )
    # Take the probability of class 1 only.
    all_probs = all_probs_[:, :, 1]

    c_score = all_probs.max(axis=1)

    # Positions of cleavage site from each model.
    possible_cleavage_sites = all_probs.argmax(axis=1) + 15
    # These positions are counted on 0 based index.
    # Make sure to add one for the 'usual' position.

    return c_score, all_probs, possible_cleavage_sites

def check_fungi(seq):
    '''
    Check if a sequence is from fungi.
    Features is the count of residues upto position 22
    '''
    seq = validate(seq)[:22]
    feat = np.array([seq.count(i) for i in 'RKNDQEHPYWSTGAMCFLVI'])
    
    classifiers = FUNGI.Classifier
    scores = np.array([clf.predict_proba([feat]) for clf in classifiers])[:, :, 1].flatten()
    return scores


def check_toxin(seq):
    '''
    Check if a sequence has toxic peptide.
    Features is hydrophobicity and SWI upto position 23
    '''
    seq = validate(seq)[:23]
    hydrop = np.array([hydrop_flex_swi[i] for i in seq])[:,0]
    swi = np.array([hydrop_flex_swi[i] for i in seq])[:,2]
    flex = np.array([hydrop_flex_swi[i] for i in seq])[:,1]
    turn = np.array([seq.count(i) for i in 'NPGS'])
    feat = np.concatenate([hydrop, swi, flex, turn])

    classifiers = TOXIN.Classifier
    scores = np.array([clf.predict_proba([feat]) for clf in classifiers])[:, :, 1].flatten()
    return scores