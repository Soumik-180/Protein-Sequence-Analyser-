def calculate_composition(seq: str) -> dict[str, float]:
    if not seq:
        return {}
    counts = {aa: seq.count(aa) for aa in set(seq)}
    total = len(seq)
    return {aa: (counts[aa] / total) * 100 for aa in counts}
