from core import (
    calculate_molecular_weight,
    calculate_composition,
    compute_charge_curve,
    calculate_pi,
    calculate_gravy,
    calculate_instability_index,
    calculate_aliphatic_index,
    calculate_hydrophobicity,
    detect_motifs,
    simple_secondary_structure_prediction
)

def format_sequence(raw_input):
    lines = raw_input.strip().splitlines()
    seq_lines = []
    found_first_header = False
    
    for line in lines:
        line = line.strip()
        if not line: continue
        if line.startswith('>'):
            if found_first_header:
                break
            found_first_header = True
            continue
        seq_lines.append(line)
        
    return ''.join(seq_lines).replace(' ', '').upper()

def run_analysis(seq):
    clean_seq = format_sequence(seq)
    if not clean_seq:
        return {"error": "Sequence is empty"}
        
    valid_aa = set("ACDEFGHIKLMNPQRSTVWY")
    if not all(aa in valid_aa for aa in clean_seq):
        return {"error": "Invalid amino acids detected in sequence. Expected 20 standard amino acids only."}
        
    return {
        "length": len(clean_seq),
        "molecular_weight": calculate_molecular_weight(clean_seq),
        "composition": calculate_composition(clean_seq),
        "charge_curve": compute_charge_curve(clean_seq),
        "isoelectric_point": calculate_pi(clean_seq),
        "gravy": calculate_gravy(clean_seq),
        "instability": calculate_instability_index(clean_seq),
        "aliphatic_index": calculate_aliphatic_index(clean_seq),
        "hydrophobicity_plot": calculate_hydrophobicity(clean_seq, window_size=9),
        "secondary_structure": simple_secondary_structure_prediction(clean_seq),
        "motifs": detect_motifs(clean_seq)
    }