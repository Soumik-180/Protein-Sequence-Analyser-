from data import AMINO_ACID_DATA


def calculate_gravy(seq: str) -> float:
    if not seq:
        return 0
    gravy_sum = sum(AMINO_ACID_DATA.get(aa, {}).get("kd", 0) for aa in seq)
    return gravy_sum / len(seq)
