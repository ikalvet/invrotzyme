import sys
import numpy as np


num2aa=[
    'ALA','ARG','ASN','ASP','CYS',
    'GLN','GLU','GLY','HIS','ILE',
    'LEU','LYS','MET','PHE','PRO',
    'SER','THR','TRP','TYR','VAL',
    ]

aa2num= {x:i for i,x in enumerate(num2aa)}

alpha_1 = list("ARNDCQEGHILKMFPSTWYV-")
aa_N_1 = {n:a for n,a in enumerate(alpha_1)}
aa_1_N = {a:n for n,a in enumerate(alpha_1)}

aa123 = {aa1: aa3 for aa1, aa3 in zip(alpha_1, num2aa)}
aa321 = {aa3: aa1 for aa1, aa3 in zip(alpha_1, num2aa)}

def N_to_AA(x):
    x = np.array(x);
    if x.ndim == 1: x = x[None]
    return ["".join([aa_N_1.get(a,"-") for a in y]) for y in x]

 
def alphabet_mapping(seq_list, alphabet_dict):
    """
    Args:
    seq_list: a list of sequences ['ABADSDAS', 'AABSDVDDV']
    Returns:
    encoded: a list of np.arrays
    """
    encoded = [[alphabet_dict[token] for token in seq] for seq in seq_list]
    return encoded


def alphabet_onehot_2_onehot(alphabet1, alphabet2):
    '''
    Args:
    alphabet1: List of amino acids in order (A characters)
    alphabet2: List of amino acids in different order
    
    Returns:
    map: AxA matrix to map one-hot encoding from alphabet1 to alphabet2
    '''
    assert len(alphabet1) == len(alphabet2), 'The alphabets must be the same length'
    
    alpha1_2_int = {aa: i for i, aa in enumerate(alphabet1)}
    alpha2_2_int = {aa: i for i, aa in enumerate(alphabet2)}
    
    A = len(alphabet1)
    map = np.zeros((A,A))
    
    for aa in alphabet1:
        j = alpha1_2_int[aa]
        i = alpha2_2_int[aa]
        map[i, j] = 1
        
    return map
        

# minimal sc atom representation (Nx8)
aa2short=[
    (" N  "," CA "," C  "," CB ",  None,  None,  None,  None), # ala
    (" N  "," CA "," C  "," CB "," CG "," CD "," NE "," CZ "), # arg
    (" N  "," CA "," C  "," CB "," CG "," OD1",  None,  None), # asn
    (" N  "," CA "," C  "," CB "," CG "," OD1",  None,  None), # asp
    (" N  "," CA "," C  "," CB "," SG ",  None,  None,  None), # cys
    (" N  "," CA "," C  "," CB "," CG "," CD "," OE1",  None), # gln
    (" N  "," CA "," C  "," CB "," CG "," CD "," OE1",  None), # glu
    (" N  "," CA "," C  ",  None,  None,  None,  None,  None), # gly
    (" N  "," CA "," C  "," CB "," CG "," ND1",  None,  None), # his
    (" N  "," CA "," C  "," CB "," CG1"," CD1",  None,  None), # ile
    (" N  "," CA "," C  "," CB "," CG "," CD1",  None,  None), # leu
    (" N  "," CA "," C  "," CB "," CG "," CD "," CE "," NZ "), # lys
    (" N  "," CA "," C  "," CB "," CG "," SD "," CE ",  None), # met
    (" N  "," CA "," C  "," CB "," CG "," CD1",  None,  None), # phe
    (" N  "," CA "," C  "," CB "," CG "," CD ",  None,  None), # pro
    (" N  "," CA "," C  "," CB "," OG ",  None,  None,  None), # ser
    (" N  "," CA "," C  "," CB "," OG1",  None,  None,  None), # thr
    (" N  "," CA "," C  "," CB "," CG "," CD1",  None,  None), # trp
    (" N  "," CA "," C  "," CB "," CG "," CD1",  None,  None), # tyr
    (" N  "," CA "," C  "," CB "," CG1",  None,  None,  None), # val
]

# full sc atom representation (Nx14)
aa2long=[
    (" N  "," CA "," C  "," O  "," CB ",  None,  None,  None,  None,  None,  None,  None,  None,  None), # ala
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD "," NE "," CZ "," NH1"," NH2",  None,  None,  None), # arg
    (" N  "," CA "," C  "," O  "," CB "," CG "," OD1"," ND2",  None,  None,  None,  None,  None,  None), # asn
    (" N  "," CA "," C  "," O  "," CB "," CG "," OD1"," OD2",  None,  None,  None,  None,  None,  None), # asp
    (" N  "," CA "," C  "," O  "," CB "," SG ",  None,  None,  None,  None,  None,  None,  None,  None), # cys
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD "," OE1"," NE2",  None,  None,  None,  None,  None), # gln
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD "," OE1"," OE2",  None,  None,  None,  None,  None), # glu
    (" N  "," CA "," C  "," O  ",  None,  None,  None,  None,  None,  None,  None,  None,  None,  None), # gly
    (" N  "," CA "," C  "," O  "," CB "," CG "," ND1"," CD2"," CE1"," NE2",  None,  None,  None,  None), # his
    (" N  "," CA "," C  "," O  "," CB "," CG1"," CG2"," CD1",  None,  None,  None,  None,  None,  None), # ile
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD1"," CD2",  None,  None,  None,  None,  None,  None), # leu
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD "," CE "," NZ ",  None,  None,  None,  None,  None), # lys
    (" N  "," CA "," C  "," O  "," CB "," CG "," SD "," CE ",  None,  None,  None,  None,  None,  None), # met
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD1"," CD2"," CE1"," CE2"," CZ ",  None,  None,  None), # phe
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD ",  None,  None,  None,  None,  None,  None,  None), # pro
    (" N  "," CA "," C  "," O  "," CB "," OG ",  None,  None,  None,  None,  None,  None,  None,  None), # ser
    (" N  "," CA "," C  "," O  "," CB "," OG1"," CG2",  None,  None,  None,  None,  None,  None,  None), # thr
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD1"," CD2"," CE2"," CE3"," NE1"," CZ2"," CZ3"," CH2"), # trp
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD1"," CD2"," CE1"," CE2"," CZ "," OH ",  None,  None), # tyr
    (" N  "," CA "," C  "," O  "," CB "," CG1"," CG2",  None,  None,  None,  None,  None,  None,  None), # val
]

# build the "alternate" sc mapping
aa2longalt=[
    (" N  "," CA "," C  "," O  "," CB ",  None,  None,  None,  None,  None,  None,  None,  None,  None), # ala
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD "," NE "," CZ "," NH2"," NH1",  None,  None,  None), # arg
    (" N  "," CA "," C  "," O  "," CB "," CG "," OD1"," ND2",  None,  None,  None,  None,  None,  None), # asn
    (" N  "," CA "," C  "," O  "," CB "," CG "," OD2"," OD1",  None,  None,  None,  None,  None,  None), # asp
    (" N  "," CA "," C  "," O  "," CB "," SG ",  None,  None,  None,  None,  None,  None,  None,  None), # cys
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD "," OE1"," NE2",  None,  None,  None,  None,  None), # gln
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD "," OE2"," OE1",  None,  None,  None,  None,  None), # glu
    (" N  "," CA "," C  "," O  ",  None,  None,  None,  None,  None,  None,  None,  None,  None,  None), # gly
    (" N  "," CA "," C  "," O  "," CB "," CG "," ND1"," CD2"," CE1"," NE2",  None,  None,  None,  None), # his
    (" N  "," CA "," C  "," O  "," CB "," CG1"," CG2"," CD1",  None,  None,  None,  None,  None,  None), # ile
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD2"," CD1",  None,  None,  None,  None,  None,  None), # leu
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD "," CE "," NZ ",  None,  None,  None,  None,  None), # lys
    (" N  "," CA "," C  "," O  "," CB "," CG "," SD "," CE ",  None,  None,  None,  None,  None,  None), # met
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD2"," CD1"," CE2"," CE1"," CZ ",  None,  None,  None), # phe
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD ",  None,  None,  None,  None,  None,  None,  None), # pro
    (" N  "," CA "," C  "," O  "," CB "," OG ",  None,  None,  None,  None,  None,  None,  None,  None), # ser
    (" N  "," CA "," C  "," O  "," CB "," OG1"," CG2",  None,  None,  None,  None,  None,  None,  None), # thr
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD1"," CD2"," CE2"," CE3"," NE1"," CZ2"," CZ3"," CH2"), # trp
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD2"," CD1"," CE2"," CE1"," CZ "," OH ",  None,  None), # tyr
    (" N  "," CA "," C  "," O  "," CB "," CG2"," CG1",  None,  None,  None,  None,  None,  None,  None), # val
]

# full sc & H atom representation (Nx22)
aa2longH = [
(' N  ', ' CA ', ' C  ', ' O  ', ' CB ', ' H  ', ' HA ', '1HB ', '2HB ', '3HB ',   None,   None,   None,   None,   None,   None,   None,   None,   None,   None,   None,   None,   None,   None), # ala
(' N  ', ' CA ', ' C  ', ' O  ', ' CB ', ' CG ', ' CD ', ' NE ', ' CZ ', ' NH1', ' NH2', ' H  ', ' HA ', '1HB ', '2HB ', '1HG ', '2HG ', '1HD ', '2HD ', ' HE ', '1HH1', '2HH1', '1HH2', '2HH2'), # arg
(' N  ', ' CA ', ' C  ', ' O  ', ' CB ', ' CG ', ' OD1', ' ND2', ' H  ', ' HA ', '1HB ', '2HB ', '1HD2', '2HD2',   None,   None,   None,   None,   None,   None,   None,   None,   None,   None), # asn
(' N  ', ' CA ', ' C  ', ' O  ', ' CB ', ' CG ', ' OD1', ' OD2', ' H  ', ' HA ', '1HB ', '2HB ',   None,   None,   None,   None,   None,   None,   None,   None,   None,   None,   None,   None), # asp
(' N  ', ' CA ', ' C  ', ' O  ', ' CB ', ' SG ', ' H  ', ' HA ', '1HB ', '2HB ', ' HG ',   None,   None,   None,   None,   None,   None,   None,   None,   None,   None,   None,   None,   None), # cys
(' N  ', ' CA ', ' C  ', ' O  ', ' CB ', ' CG ', ' CD ', ' OE1', ' NE2', ' H  ', ' HA ', '1HB ', '2HB ', '1HG ', '2HG ', '1HE2', '2HE2',   None,   None,   None,   None,   None,   None,   None), # gln
(' N  ', ' CA ', ' C  ', ' O  ', ' CB ', ' CG ', ' CD ', ' OE1', ' OE2', ' H  ', ' HA ', '1HB ', '2HB ', '1HG ', '2HG ',   None,   None,   None,   None,   None,   None,   None,   None,   None), # glu
(' N  ', ' CA ', ' C  ', ' O  ', ' H  ', '1HA ', '2HA ',   None,   None,   None,   None,   None,   None,   None,   None,   None,   None,   None,   None,   None,   None,   None,   None,   None), # gly
(' N  ', ' CA ', ' C  ', ' O  ', ' CB ', ' CG ', ' ND1', ' CD2', ' CE1', ' NE2', ' H  ', ' HA ', '1HB ', '2HB ', ' HD2', ' HE1', ' HE2',   None,   None,   None,   None,   None,   None,   None), # his
(' N  ', ' CA ', ' C  ', ' O  ', ' CB ', ' CG1', ' CG2', ' CD1', ' H  ', ' HA ', ' HB ', '1HG1', '2HG1', '1HG2', '2HG2', '3HG2', '1HD1', '2HD1', '3HD1',   None,   None,   None,   None,   None), # ile
(' N  ', ' CA ', ' C  ', ' O  ', ' CB ', ' CG ', ' CD1', ' CD2', ' H  ', ' HA ', '1HB ', '2HB ', ' HG ', '1HD1', '2HD1', '3HD1', '1HD2', '2HD2', '3HD2',   None,   None,   None,   None,   None), # leu
(' N  ', ' CA ', ' C  ', ' O  ', ' CB ', ' CG ', ' CD ', ' CE ', ' NZ ', ' H  ', ' HA ', '1HB ', '2HB ', '1HG ', '2HG ', '1HD ', '2HD ', '1HE ', '2HE ', '1HZ ', '2HZ ', '3HZ ',   None,   None), # lys
(' N  ', ' CA ', ' C  ', ' O  ', ' CB ', ' CG ', ' SD ', ' CE ', ' H  ', ' HA ', '1HB ', '2HB ', '1HG ', '2HG ', '1HE ', '2HE ', '3HE ',   None,   None,   None,   None,   None,   None,   None), # met
(' N  ', ' CA ', ' C  ', ' O  ', ' CB ', ' CG ', ' CD1', ' CD2', ' CE1', ' CE2', ' CZ ', ' H  ', ' HA ', '1HB ', '2HB ', ' HD1', ' HD2', ' HE1', ' HE2', ' HZ ',   None,   None,   None,   None), # phe
(' N  ', ' CA ', ' C  ', ' O  ', ' CB ', ' CG ', ' CD ', ' NV ', ' HA ', '1HB ', '2HB ', '1HG ', '2HG ', '1HD ', '2HD ',   None,   None,   None,   None,   None,   None,   None,   None,   None), # pro
(' N  ', ' CA ', ' C  ', ' O  ', ' CB ', ' OG ', ' H  ', ' HA ', '1HB ', '2HB ', ' HG ',   None,   None,   None,   None,   None,   None,   None,   None,   None,   None,   None,   None,   None), # ser
(' N  ', ' CA ', ' C  ', ' O  ', ' CB ', ' OG1', ' CG2', ' H  ', ' HA ', ' HB ', ' HG1', '1HG2', '2HG2', '3HG2',   None,   None,   None,   None,   None,   None,   None,   None,   None,   None), # thr
(' N  ', ' CA ', ' C  ', ' O  ', ' CB ', ' CG ', ' CD1', ' CD2', ' NE1', ' CE2', ' CE3', ' CZ2', ' CZ3', ' CH2', ' H  ', ' HA ', '1HB ', '2HB ', ' HD1', ' HE1', ' HE3', ' HZ2',  'HZ3',  'HH2'), # trp
(' N  ', ' CA ', ' C  ', ' O  ', ' CB ', ' CG ', ' CD1', ' CD2', ' CE1', ' CE2', ' CZ ', ' OH ', ' H  ', ' HA ', '1HB ', '2HB ', ' HD1', ' HD2', ' HE1', ' HE2', ' HH ',   None,   None,   None), # tyr
(' N  ', ' CA ', ' C  ', ' O  ', ' CB ', ' CG1', ' CG2', ' H  ', ' HA ', ' HB ', '1HG1', '2HG1', '3HG1', '1HG2', '2HG2', '3HG2',   None,   None,   None,   None,   None,   None,   None,   None)  # val
]

aa2longH_Nterm = [
(' N  ', ' CA ', ' C  ', ' O  ', ' CB ', '1H  ', '2H  ', '3H  ', ' HA ', '1HB ', '2HB ', '3HB ',   None,   None,   None,   None,   None,   None,   None,   None,   None,   None,   None,   None), # ala
(' N  ', ' CA ', ' C  ', ' O  ', ' CB ', ' CG ', ' CD ', ' NE ', ' CZ ', ' NH1', ' NH2', '1H  ', '2H  ', '3H  ', ' HA ', '1HB ', '2HB ', '1HG ', '2HG ', '1HD ', '2HD ', ' HE ', '1HH1', '2HH1', '1HH2', '2HH2'), # arg
(' N  ', ' CA ', ' C  ', ' O  ', ' CB ', ' CG ', ' OD1', ' ND2', '1H  ', '2H  ', '3H  ', ' HA ', '1HB ', '2HB ', '1HD2', '2HD2',   None,   None,   None,   None,   None,   None,   None,   None,   None,   None), # asn
(' N  ', ' CA ', ' C  ', ' O  ', ' CB ', ' CG ', ' OD1', ' OD2', '1H  ', '2H  ', '3H  ', ' HA ', '1HB ', '2HB ',   None,   None,   None,   None,   None,   None,   None,   None,   None,   None,   None,   None), # asp
(' N  ', ' CA ', ' C  ', ' O  ', ' CB ', ' SG ', '1H  ', '2H  ', '3H  ', ' HA ', '1HB ', '2HB ', ' HG ',   None,   None,   None,   None,   None,   None,   None,   None,   None,   None,   None,   None,   None), # cys
(' N  ', ' CA ', ' C  ', ' O  ', ' CB ', ' CG ', ' CD ', ' OE1', ' NE2', '1H  ', '2H  ', '3H  ', ' HA ', '1HB ', '2HB ', '1HG ', '2HG ', '1HE2', '2HE2',   None,   None,   None,   None,   None,   None,   None), # gln
(' N  ', ' CA ', ' C  ', ' O  ', ' CB ', ' CG ', ' CD ', ' OE1', ' OE2', '1H  ', '2H  ', '3H  ', ' HA ', '1HB ', '2HB ', '1HG ', '2HG ',   None,   None,   None,   None,   None,   None,   None,   None,   None), # glu
(' N  ', ' CA ', ' C  ', ' O  ', '1H  ', '2H  ', '3H  ', '1HA ', '2HA ',   None,   None,   None,   None,   None,   None,   None,   None,   None,   None,   None,   None,   None,   None,   None,   None,   None), # gly
(' N  ', ' CA ', ' C  ', ' O  ', ' CB ', ' CG ', ' ND1', ' CD2', ' CE1', ' NE2', '1H  ', '2H  ', '3H  ', ' HA ', '1HB ', '2HB ', ' HD2', ' HE1', ' HE2',   None,   None,   None,   None,   None,   None,   None), # his
(' N  ', ' CA ', ' C  ', ' O  ', ' CB ', ' CG1', ' CG2', ' CD1', '1H  ', '2H  ', '3H  ', ' HA ', ' HB ', '1HG1', '2HG1', '1HG2', '2HG2', '3HG2', '1HD1', '2HD1', '3HD1',   None,   None,   None,   None,   None), # ile
(' N  ', ' CA ', ' C  ', ' O  ', ' CB ', ' CG ', ' CD1', ' CD2', '1H  ', '2H  ', '3H  ', ' HA ', '1HB ', '2HB ', ' HG ', '1HD1', '2HD1', '3HD1', '1HD2', '2HD2', '3HD2',   None,   None,   None,   None,   None), # leu
(' N  ', ' CA ', ' C  ', ' O  ', ' CB ', ' CG ', ' CD ', ' CE ', ' NZ ', '1H  ', '2H  ', '3H  ', ' HA ', '1HB ', '2HB ', '1HG ', '2HG ', '1HD ', '2HD ', '1HE ', '2HE ', '1HZ ', '2HZ ', '3HZ ',   None,   None), # lys
(' N  ', ' CA ', ' C  ', ' O  ', ' CB ', ' CG ', ' SD ', ' CE ', '1H  ', '2H  ', '3H  ', ' HA ', '1HB ', '2HB ', '1HG ', '2HG ', '1HE ', '2HE ', '3HE ',   None,   None,   None,   None,   None,   None,   None), # met
(' N  ', ' CA ', ' C  ', ' O  ', ' CB ', ' CG ', ' CD1', ' CD2', ' CE1', ' CE2', ' CZ ', '1H  ', '2H  ', '3H  ', ' HA ', '1HB ', '2HB ', ' HD1', ' HD2', ' HE1', ' HE2', ' HZ ',   None,   None,   None,   None), # phe
(' N  ', ' CA ', ' C  ', ' O  ', ' CB ', ' CG ', ' CD ', ' NV ', 'CAV ', ' HA ', '1HB ', '2HB ', '1HG ', '2HG ', '1HD ', '2HD ', '1H  ', '2H  ',   None,   None,   None,   None,   None,   None,   None,   None), # pro
(' N  ', ' CA ', ' C  ', ' O  ', ' CB ', ' OG ', '1H  ', '2H  ', '3H  ', ' HA ', '1HB ', '2HB ', ' HG ',   None,   None,   None,   None,   None,   None,   None,   None,   None,   None,   None,   None,   None), # ser
(' N  ', ' CA ', ' C  ', ' O  ', ' CB ', ' OG1', ' CG2', '1H  ', '2H  ', '3H  ', ' HA ', ' HB ', ' HG1', '1HG2', '2HG2', '3HG2',   None,   None,   None,   None,   None,   None,   None,   None,   None,   None), # thr
(' N  ', ' CA ', ' C  ', ' O  ', ' CB ', ' CG ', ' CD1', ' CD2', ' NE1', ' CE2', ' CE3', ' CZ2', ' CZ3', ' CH2', '1H  ', '2H  ', '3H  ', ' HA ', '1HB ', '2HB ', ' HD1', ' HE1', ' HE3', ' HZ2',  'HZ3',  'HH2'), # trp
(' N  ', ' CA ', ' C  ', ' O  ', ' CB ', ' CG ', ' CD1', ' CD2', ' CE1', ' CE2', ' CZ ', ' OH ', '1H  ', '2H  ', '3H  ', ' HA ', '1HB ', '2HB ', ' HD1', ' HD2', ' HE1', ' HE2', ' HH ',   None,   None,   None), # tyr
(' N  ', ' CA ', ' C  ', ' O  ', ' CB ', ' CG1', ' CG2', '1H  ', '2H  ', '3H  ', ' HA ', ' HB ', '1HG1', '2HG1', '3HG1', '1HG2', '2HG2', '3HG2',   None,   None,   None,   None,   None,   None,   None,   None)  # val
]

aa2longH_Cterm = [
(' N  ', ' CA ', ' C  ', ' O  ', ' OXT', ' CB ', ' H  ', ' HA ', '1HB ', '2HB ', '3HB ',   None,   None,   None,   None,   None,   None,   None,   None,   None,   None,   None,   None,   None,   None), # ala
(' N  ', ' CA ', ' C  ', ' O  ', ' OXT', ' CB ', ' CG ', ' CD ', ' NE ', ' CZ ', ' NH1', ' NH2', ' H  ', ' HA ', '1HB ', '2HB ', '1HG ', '2HG ', '1HD ', '2HD ', ' HE ', '1HH1', '2HH1', '1HH2', '2HH2'), # arg
(' N  ', ' CA ', ' C  ', ' O  ', ' OXT', ' CB ', ' CG ', ' OD1', ' ND2', ' H  ', ' HA ', '1HB ', '2HB ', '1HD2', '2HD2',   None,   None,   None,   None,   None,   None,   None,   None,   None,   None), # asn
(' N  ', ' CA ', ' C  ', ' O  ', ' OXT', ' CB ', ' CG ', ' OD1', ' OD2', ' H  ', ' HA ', '1HB ', '2HB ',   None,   None,   None,   None,   None,   None,   None,   None,   None,   None,   None,   None), # asp
(' N  ', ' CA ', ' C  ', ' O  ', ' OXT', ' CB ', ' SG ', ' H  ', ' HA ', '1HB ', '2HB ', ' HG ',   None,   None,   None,   None,   None,   None,   None,   None,   None,   None,   None,   None,   None), # cys
(' N  ', ' CA ', ' C  ', ' O  ', ' OXT', ' CB ', ' CG ', ' CD ', ' OE1', ' NE2', ' H  ', ' HA ', '1HB ', '2HB ', '1HG ', '2HG ', '1HE2', '2HE2',   None,   None,   None,   None,   None,   None,   None), # gln
(' N  ', ' CA ', ' C  ', ' O  ', ' OXT', ' CB ', ' CG ', ' CD ', ' OE1', ' OE2', ' H  ', ' HA ', '1HB ', '2HB ', '1HG ', '2HG ',   None,   None,   None,   None,   None,   None,   None,   None,   None), # glu
(' N  ', ' CA ', ' C  ', ' O  ', ' OXT', ' H  ', '1HA ', '2HA ',   None,   None,   None,   None,   None,   None,   None,   None,   None,   None,   None,   None,   None,   None,   None,   None,   None), # gly
(' N  ', ' CA ', ' C  ', ' O  ', ' OXT', ' CB ', ' CG ', ' ND1', ' CD2', ' CE1', ' NE2', ' H  ', ' HA ', '1HB ', '2HB ', ' HD2', ' HE1', ' HE2',   None,   None,   None,   None,   None,   None,   None), # his
(' N  ', ' CA ', ' C  ', ' O  ', ' OXT', ' CB ', ' CG1', ' CG2', ' CD1', ' H  ', ' HA ', ' HB ', '1HG1', '2HG1', '1HG2', '2HG2', '3HG2', '1HD1', '2HD1', '3HD1',   None,   None,   None,   None,   None), # ile
(' N  ', ' CA ', ' C  ', ' O  ', ' OXT', ' CB ', ' CG ', ' CD1', ' CD2', ' H  ', ' HA ', '1HB ', '2HB ', ' HG ', '1HD1', '2HD1', '3HD1', '1HD2', '2HD2', '3HD2',   None,   None,   None,   None,   None), # leu
(' N  ', ' CA ', ' C  ', ' O  ', ' OXT', ' CB ', ' CG ', ' CD ', ' CE ', ' NZ ', ' H  ', ' HA ', '1HB ', '2HB ', '1HG ', '2HG ', '1HD ', '2HD ', '1HE ', '2HE ', '1HZ ', '2HZ ', '3HZ ',   None,   None), # lys
(' N  ', ' CA ', ' C  ', ' O  ', ' OXT', ' CB ', ' CG ', ' SD ', ' CE ', ' H  ', ' HA ', '1HB ', '2HB ', '1HG ', '2HG ', '1HE ', '2HE ', '3HE ',   None,   None,   None,   None,   None,   None,   None), # met
(' N  ', ' CA ', ' C  ', ' O  ', ' OXT', ' CB ', ' CG ', ' CD1', ' CD2', ' CE1', ' CE2', ' CZ ', ' H  ', ' HA ', '1HB ', '2HB ', ' HD1', ' HD2', ' HE1', ' HE2', ' HZ ',   None,   None,   None,   None), # phe
(' N  ', ' CA ', ' C  ', ' O  ', ' OXT', ' CB ', ' CG ', ' CD ', ' NV ', ' HA ', '1HB ', '2HB ', '1HG ', '2HG ', '1HD ', '2HD ',   None,   None,   None,   None,   None,   None,   None,   None,   None), # pro
(' N  ', ' CA ', ' C  ', ' O  ', ' OXT', ' CB ', ' OG ', ' H  ', ' HA ', '1HB ', '2HB ', ' HG ',   None,   None,   None,   None,   None,   None,   None,   None,   None,   None,   None,   None,   None), # ser
(' N  ', ' CA ', ' C  ', ' O  ', ' OXT', ' CB ', ' OG1', ' CG2', ' H  ', ' HA ', ' HB ', ' HG1', '1HG2', '2HG2', '3HG2',   None,   None,   None,   None,   None,   None,   None,   None,   None,   None), # thr
(' N  ', ' CA ', ' C  ', ' O  ', ' OXT', ' CB ', ' CG ', ' CD1', ' CD2', ' NE1', ' CE2', ' CE3', ' CZ2', ' CZ3', ' CH2', ' H  ', ' HA ', '1HB ', '2HB ', ' HD1', ' HE1', ' HE3', ' HZ2',  'HZ3',  'HH2'), # trp
(' N  ', ' CA ', ' C  ', ' O  ', ' OXT', ' CB ', ' CG ', ' CD1', ' CD2', ' CE1', ' CE2', ' CZ ', ' OH ', ' H  ', ' HA ', '1HB ', '2HB ', ' HD1', ' HD2', ' HE1', ' HE2', ' HH ',   None,   None,   None), # tyr
(' N  ', ' CA ', ' C  ', ' O  ', ' OXT', ' CB ', ' CG1', ' CG2', ' H  ', ' HA ', ' HB ', '1HG1', '2HG1', '3HG1', '1HG2', '2HG2', '3HG2',   None,   None,   None,   None,   None,   None,   None,   None)  # val
]


# build "deterministic" atoms
# see notebook (se3_experiments.ipynb for derivation)
aa2frames=[
    [], # ala
    [   # arg
        [' NH1', ' CZ ', ' NE ', ' CD ', [-0.7218378782272339, 1.0856682062149048, -0.006118079647421837]],
        [' NH2', ' CZ ', ' NE ', ' CD ', [-0.6158039569854736, -1.1400136947631836, 0.006467342376708984]]],
    [   # asn
        [' ND2', ' CG ', ' CB ', ' OD1', [-0.6304131746292114, -1.1431225538253784, 0.02364802360534668]]],
    [   # asp
        [' OD2', ' CG ', ' CB ', ' OD1', [-0.5972501039505005, -1.0955055952072144, 0.04530305415391922]]],
    [], # cys
    [   # gln
        [' NE2', ' CD ', ' CG ', ' OE1', [-0.6558755040168762, -1.1324536800384521, 0.026521772146224976]]],
    [   # glu
        [' OE2', ' CD ', ' CG ', ' OE1', [-0.5578438639640808, -1.1161314249038696, -0.015464287251234055]]],
    [], # gly
    [   # his
        [' CD2', ' CG ', ' CB ', ' ND1', [-0.7502505779266357, -1.1680538654327393, 0.0005368441343307495]],
        [' CE1', ' CG ', ' CB ', ' ND1', [-2.0262467861175537, 0.539483368396759, -0.004495501518249512]],
        [' NE2', ' CG ', ' CB ', ' ND1', [-2.0761325359344482, -0.8199722766876221, -0.0018703639507293701]]],
    [   # ile
        [' CG2', ' CB ', ' CA ', ' CG1', [-0.6059935688972473, -0.8108057379722595, 1.1861376762390137]]],
    [   # leu
        [' CD2', ' CG ', ' CB ', ' CD1', [-0.5942193269729614, -0.7693282961845398, -1.1914138793945312]]],
    [], # lys
    [], # met
    [   # phe
        [' CD2', ' CG ', ' CB ', ' CD1', [-0.7164441347122192, -1.197853446006775, 0.06416648626327515]],
        [' CE1', ' CG ', ' CB ', ' CD1', [-2.0785865783691406, 1.2366485595703125, 0.08100450038909912]],
        [' CE2', ' CG ', ' CB ', ' CD1', [-2.107091188430786, -1.178497076034546, 0.13524535298347473]],
        [' CZ ', ' CG ', ' CB ', ' CD1', [-2.786630630493164, 0.03873880207538605, 0.14633776247501373]]],
    [], # pro
    [], # ser
    [   # thr
        [' CG2', ' CB ', ' CA ', ' OG1', [-0.6842088103294373, -0.6709619164466858, 1.2105456590652466]]],
    [   # trp
        [' CD2', ' CG ', ' CB ', ' CD1', [-0.8550368547439575, -1.0790592432022095, 0.09017711877822876]],
        [' NE1', ' CG ', ' CB ', ' CD1', [-2.1863200664520264, 0.8064242601394653, 0.08350661396980286]],
        [' CE2', ' CG ', ' CB ', ' CD1', [-2.1801204681396484, -0.5795643329620361, 0.14015203714370728]],
        [' CE3', ' CG ', ' CB ', ' CD1', [-0.605582594871521, -2.4733362197875977, 0.16200461983680725]],
        [' CE2', ' CG ', ' CB ', ' CD1', [-2.1801204681396484, -0.5795643329620361, 0.14015203714370728]],
        [' CZ2', ' CG ', ' CB ', ' CD1', [-3.2672977447509766, -1.473116159439087, 0.250858873128891]],
        [' CZ3', ' CG ', ' CB ', ' CD1', [-1.6969941854476929, -3.3360071182250977, 0.264143705368042]],
        [' CH2', ' CG ', ' CB ', ' CD1', [-3.009331703186035, -2.8451972007751465, 0.3059283494949341]]],
    [   # tyr
        [' CD2', ' CG ', ' CB ', ' CD1', [-0.69439297914505, -1.2123756408691406, -0.009198814630508423]],
        [' CE1', ' CG ', ' CB ', ' CD1', [-2.104464054107666, 1.1910505294799805, -0.014679580926895142]],
        [' CE2', ' CG ', ' CB ', ' CD1', [-2.0857787132263184, -1.2231677770614624, -0.024517983198165894]],
        [' CZ ', ' CG ', ' CB ', ' CD1', [-2.7897322177886963, -0.021470561623573303, -0.026979409158229828]],
        [' OH ', ' CG ', ' CB ', ' CD1', [-4.1559271812438965, -0.029129385948181152, -0.044720835983753204]]],
    [   # val
        [' CG2', ' CB ', ' CA ', ' CG1', [-0.6258467435836792, -0.7654698491096497, -1.1894742250442505]]],
]

# O from frame (C,N-1,CA)
bb2oframe=[-0.5992066264152527, -1.0820008516311646, 0.0001476481556892395]

# build the mapping from indices in reduced representation to 
# indices in the full representation
#  N x 14 x 6 = <base-idx | parent-idx | gparent-idx | x | y | z >
#    base-idx < 0 ==> no atom
#    xyz = 0 ==> no mapping
short2long = np.zeros((20,14,6))
for i in range(20):
    i_s, i_l = aa2short[i],aa2long[i]
    for j,a in enumerate(i_l):
        # case 1: if no atom defined, blank
        if (a is None):
            short2long[i,j,0] = -1
        # case 2: atom is a base atom
        elif (a in i_s):
            short2long[i,j,0] = i_s.index(a)
            if (short2long[i,j,0] == 0):
                short2long[i,j,1] = 1
                short2long[i,j,2] = 2
            else:
                short2long[i,j,1] = 0
                if (short2long[i,j,0] == 1):
                    short2long[i,j,2] = 2
                else:
                    short2long[i,j,2] = 1
        # case 3: atom is ' O  '
        elif (a == " O  "):
            short2long[i,j,0] = 2
            short2long[i,j,1] = 0 #Nprev (will pre-roll N as nothing else needs it)
            short2long[i,j,2] = 1
            short2long[i,j,3:] = np.array(bb2oframe)
        # case 4: build this atom
        else:
            i_f = aa2frames[i]
            names = [f[0] for f in i_f]
            idx = names.index(a)
            short2long[i,j,0] = i_s.index(i_f[idx][1])
            short2long[i,j,1] = i_s.index(i_f[idx][2])
            short2long[i,j,2] = i_s.index(i_f[idx][3])
            short2long[i,j,3:] = np.array(i_f[idx][4])

# build the mapping from atoms in the full rep (Nx14) to the "alternate" rep
long2alt = np.zeros((20,14))
for i in range(20):
    i_l, i_lalt = aa2long[i],  aa2longalt[i]
    for j,a in enumerate(i_l):
        if (a is None):
            long2alt[i,j] = j
        else:
            long2alt[i,j] = i_lalt.index(a)

