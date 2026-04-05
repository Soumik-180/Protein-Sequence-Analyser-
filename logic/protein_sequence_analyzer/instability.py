from typing import Union

from data import DIWV


def calculate_instability_index(seq: str) -> dict[str, Union[float, bool]]:
    """Calculates instability using the full DIWV matrix."""
    if len(seq) < 2:
        return {"score": 0, "stable": True}

    score = 0.0
    for i in range(len(seq) - 1):
        pair_start = seq[i]
        pair_end = seq[i + 1]
        score += DIWV.get(pair_start, {}).get(pair_end, 0.0)

    final_score = (10.0 / len(seq)) * score
    return {"score": final_score, "stable": final_score < 40}
