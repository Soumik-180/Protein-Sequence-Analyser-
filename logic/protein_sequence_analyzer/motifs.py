import re
from typing import Union


def detect_motifs(seq: str) -> list[dict[str, Union[int, str]]]:
    motifs: list[dict[str, Union[int, str]]] = []

    for m in re.finditer(r"(?=N[^P][ST][^P])", seq):
        motifs.append({"type": "N-glycosylation", "position": m.start() + 1})

    for m in re.finditer(r"(?=[ST]..[DE])", seq):
        motifs.append({"type": "Casein kinase II phosphorylation", "position": m.start() + 1})

    return motifs
