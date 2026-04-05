from logic.protein_sequence_analyzer.molecular_weight import calculate_molecular_weight
from logic.protein_sequence_analyzer.composition import calculate_composition
from logic.protein_sequence_analyzer.charge_curve import compute_charge_curve
from logic.protein_sequence_analyzer.isoelectric_point import calculate_pi
from logic.protein_sequence_analyzer.gravy import calculate_gravy
from logic.protein_sequence_analyzer.instability import calculate_instability_index
from logic.protein_sequence_analyzer.aliphatic_index import calculate_aliphatic_index
from logic.protein_sequence_analyzer.hydrophobicity import calculate_hydrophobicity
from logic.protein_sequence_analyzer.motifs import detect_motifs
from logic.protein_sequence_analyzer.secondary_structure import simple_secondary_structure_prediction

from logic.amino_acid_explorer.normalize import normalize_amino_acid_query
from logic.amino_acid_explorer.report import build_amino_acid_report

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


def run_amino_acid_explorer(query: str):
    code, error = normalize_amino_acid_query(query)
    if error:
        return {"error": error}

    try:
        return build_amino_acid_report(code)
    except Exception:
        return {"error": "Failed to build amino acid report."}
