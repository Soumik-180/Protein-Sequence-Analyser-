from data import AMINO_ACID_DATA


def calculate_hydrophobicity(seq: str, window_size: int = 9) -> list[dict[str, float]]:
    if len(seq) < window_size:
        return []

    plot_data: list[dict[str, float]] = []
    for i in range(len(seq) - window_size + 1):
        window = seq[i : i + window_size]
        score = sum(AMINO_ACID_DATA.get(aa, {}).get("kd", 0) for aa in window) / window_size
        plot_data.append({"position": i + window_size // 2 + 1, "score": score})

    return plot_data
