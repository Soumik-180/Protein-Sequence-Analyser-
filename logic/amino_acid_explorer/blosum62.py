from data import BLOSUM62


def get_blosum62_row(code: str) -> dict[str, int]:
    row = BLOSUM62.get(code, {})
    return dict(row)


def get_blosum62_score(base: str, other: str):
    return BLOSUM62.get(base, {}).get(other)
