def simple_secondary_structure_prediction(seq: str) -> list[str]:
    structure: list[str] = []
    for aa in seq:
        if aa in ["E", "A", "L", "M", "Q"]:
            structure.append("Alpha Helix")
        elif aa in ["V", "I", "Y", "C", "W", "F"]:
            structure.append("Beta Sheet")
        else:
            structure.append("Coil")
    return structure
