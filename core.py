import re
import numpy as np
from data import AMINO_ACID_DATA, BLOSUM62, DIWV

def calculate_molecular_weight(seq):
    weight = sum(AMINO_ACID_DATA.get(aa, {}).get('mw', 0) for aa in seq)
    if len(seq) > 1:
        weight -= (len(seq) - 1) * 18.015
    return weight

def calculate_composition(seq):
    if not seq: return {}
    counts = {aa: seq.count(aa) for aa in set(seq)}
    total = len(seq)
    return {aa: (counts[aa] / total) * 100 for aa in counts}

PKA_VALUES = {
    'N_terminus': 9.60, 'C_terminus': 3.60,
    'D': 3.87, 'E': 4.25, 'C': 8.36, 
    'H': 6.04, 'K': 10.45, 'R': 12.48, 'Y': 10.46
}

def _charge_at_ph(seq, pH):
    charge = 1 / (1 + 10 ** (pH - PKA_VALUES['N_terminus']))
    charge -= 1 / (1 + 10 ** (PKA_VALUES['C_terminus'] - pH))
    for aa in seq:
        if aa in ['D', 'E', 'C', 'Y']:
            charge -= 1 / (1 + 10 ** (PKA_VALUES[aa] - pH))
        elif aa in ['K', 'R', 'H']:
            charge += 1 / (1 + 10 ** (pH - PKA_VALUES[aa]))
    return charge

def compute_charge_curve(seq):
    ph_values = np.arange(0.0, 14.5, 0.5)
    charges = [_charge_at_ph(seq, ph) for ph in ph_values]
    return {"ph": ph_values.tolist(), "charge": charges}

def calculate_pi(seq):
    low_ph, high_ph, pi = 0.0, 14.0, 7.0
    for _ in range(100):
        pi = (low_ph + high_ph) / 2
        if _charge_at_ph(seq, pi) > 0:
            low_ph = pi
        else:
            high_ph = pi
    return pi

def calculate_gravy(seq):
    if not seq:
        return 0
    # Consistently uses 'kd'
    gravy_sum = sum(AMINO_ACID_DATA.get(aa, {}).get('kd', 0) for aa in seq)
    return gravy_sum / len(seq)

def calculate_instability_index(seq):
    """Calculates instability using the full DIWV matrix"""
    if len(seq) < 2:
        return {'score': 0, 'stable': True}
    
    score = 0.0
    for i in range(len(seq) - 1):
        pair_start = seq[i]
        pair_end = seq[i+1]
        score += DIWV.get(pair_start, {}).get(pair_end, 0.0)
        
    final_score = (10.0 / len(seq)) * score
    return {'score': final_score, 'stable': final_score < 40}

def calculate_aliphatic_index(seq):
    if not seq:
        return 0
    # Guard against unknown codes
    valid_seq = [aa for aa in seq if aa in set("ACDEFGHIKLMNPQRSTVWY")]
    if not valid_seq: return 0
    
    a_count = valid_seq.count('A')
    v_count = valid_seq.count('V')
    i_count = valid_seq.count('I')
    l_count = valid_seq.count('L')
    return ((a_count + 2.9 * v_count + 3.9 * (i_count + l_count)) / len(valid_seq)) * 100

def calculate_hydrophobicity(seq, window_size=9):
    if len(seq) < window_size:
        return []
    plot_data = []
    for i in range(len(seq) - window_size + 1):
        window = seq[i:i + window_size]
        score = sum(AMINO_ACID_DATA.get(aa, {}).get('kd', 0) for aa in window) / window_size
        plot_data.append({'position': i + window_size // 2 + 1, 'score': score})
    return plot_data

def detect_motifs(seq):
    motifs = []
    for m in re.finditer(r'(?=N[^P][ST][^P])', seq):
        motifs.append({"type": "N-glycosylation", "position": m.start() + 1})
    for m in re.finditer(r'(?=[ST]..[DE])', seq):
        motifs.append({"type": "Casein kinase II phosphorylation", "position": m.start() + 1})
    return motifs

def simple_secondary_structure_prediction(seq):
    structure = []
    for aa in seq:
        if aa in ['E', 'A', 'L', 'M', 'Q']:
            structure.append('Alpha Helix')
        elif aa in ['V', 'I', 'Y', 'C', 'W', 'F']:
            structure.append('Beta Sheet')
        else:
            structure.append('Coil')
    return structure

def smith_waterman(seq1, seq2, gap_open=-11, gap_ext=-1):
    n, m = len(seq1), len(seq2)
    if n == 0 or m == 0:
        return {'score': 0, 'similarity': 0, 'matches': 0, 'alignLen': 0, 'aligned_seq1': '', 'aligned_seq2': ''}

    H = [[0] * (m + 1) for _ in range(n + 1)]
    E = [[0] * (m + 1) for _ in range(n + 1)]
    F = [[0] * (m + 1) for _ in range(n + 1)]
    
    # State flags: 0=Stop, 1=H(diag), 2=E(up), 3=F(left)
    TB_H = [[0] * (m + 1) for _ in range(n + 1)]
    TB_E = [[0] * (m + 1) for _ in range(n + 1)]
    TB_F = [[0] * (m + 1) for _ in range(n + 1)]

    max_score = 0
    max_i, max_j = 0, 0

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            match_score = BLOSUM62.get(seq1[i-1], {}).get(seq2[j-1], -4)
            
            # E (gap in seq2 / vertical)
            e_open = H[i-1][j] + gap_open
            e_ext = E[i-1][j] + gap_ext
            if e_open >= e_ext:
                E[i][j] = e_open
                TB_E[i][j] = 1 # from H
            else:
                E[i][j] = e_ext
                TB_E[i][j] = 2 # from E

            # F (gap in seq1 / horizontal)
            f_open = H[i][j-1] + gap_open
            f_ext = F[i][j-1] + gap_ext
            if f_open >= f_ext:
                F[i][j] = f_open
                TB_F[i][j] = 1 # from H
            else:
                F[i][j] = f_ext
                TB_F[i][j] = 3 # from F

            # H
            diag = H[i-1][j-1] + match_score
            scores = (0, diag, E[i][j], F[i][j])
            best = max(scores)
            H[i][j] = best
            
            if best == 0:
                TB_H[i][j] = 0
            elif best == diag:
                TB_H[i][j] = 1
            elif best == E[i][j]:
                TB_H[i][j] = 2
            elif best == F[i][j]:
                TB_H[i][j] = 3

            if H[i][j] > max_score:
                max_score = H[i][j]
                max_i, max_j = i, j

    i, j = max_i, max_j
    matches, align_len = 0, 0
    al1, al2 = [], []
    state = 1 # Start traceback in H state
    
    while i > 0 and j > 0:
        if state == 1: # In H matrix
            if TB_H[i][j] == 0:
                break
            tr = TB_H[i][j]
            if tr == 1:
                align_len += 1
                if seq1[i-1] == seq2[j-1]: matches += 1
                al1.append(seq1[i-1])
                al2.append(seq2[j-1])
                i -= 1; j -= 1
            elif tr == 2:
                state = 2
            elif tr == 3:
                state = 3
        elif state == 2: # In E matrix
            align_len += 1
            al1.append(seq1[i-1])
            al2.append('-')
            state = TB_E[i][j]
            i -= 1
        elif state == 3: # In F matrix
            align_len += 1
            al1.append('-')
            al2.append(seq2[j-1])
            state = TB_F[i][j]
            j -= 1

    similarity = (matches / align_len * 100) if align_len > 0 else 0
    return {
        'score': max_score, 
        'similarity': similarity, 
        'matches': matches, 
        'alignLen': align_len,
        'aligned_seq1': ''.join(reversed(al1)),
        'aligned_seq2': ''.join(reversed(al2))
    }
