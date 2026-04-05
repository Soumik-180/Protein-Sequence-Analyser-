from data import AMINO_ACID_DATA


def normalize_amino_acid_query(query: str):
    """Normalize an amino-acid query to a 1-letter code.

    Accepts:
    - 1-letter code (e.g., "A")
    - 3-letter code (e.g., "Ala")
    - Full name (e.g., "Alanine")

    Returns:
        (code, error) where exactly one of them is non-None.
    """
    if query is None:
        return None, "Please enter an amino acid."

    raw = str(query).strip()
    if not raw:
        return None, "Please enter an amino acid."

    q_upper = raw.upper()
    if len(q_upper) == 1 and q_upper in AMINO_ACID_DATA:
        return q_upper, None

    three_to_one = {
        (props.get("three") or "").upper(): one
        for one, props in AMINO_ACID_DATA.items()
        if props.get("three")
    }
    if q_upper in three_to_one:
        return three_to_one[q_upper], None

    q_lower = raw.lower()
    name_to_one = {
        (props.get("name") or "").lower(): one
        for one, props in AMINO_ACID_DATA.items()
        if props.get("name")
    }
    if q_lower in name_to_one:
        return name_to_one[q_lower], None

    # Convenience: allow "aspartic acid" / "glutamic acid"
    if q_lower.endswith("ic acid"):
        alt = q_lower.replace("ic acid", "ate")
        if alt in name_to_one:
            return name_to_one[alt], None

    return (
        None,
        "Unknown amino acid. Use 1-letter (A), 3-letter (Ala), or full name (Alanine).",
    )
