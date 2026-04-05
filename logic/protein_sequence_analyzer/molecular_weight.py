from data import AMINO_ACID_DATA


def calculate_molecular_weight(seq: str) -> float:
    weight = sum(AMINO_ACID_DATA.get(aa, {}).get("mw", 0) for aa in seq)
    if len(seq) > 1:
        # peptide bond formation releases water
        weight -= (len(seq) - 1) * 18.015
    return weight
