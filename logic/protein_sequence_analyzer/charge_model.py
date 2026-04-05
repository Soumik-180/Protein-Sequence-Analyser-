PKA_VALUES = {
    "N_terminus": 9.60,
    "C_terminus": 3.60,
    "D": 3.87,
    "E": 4.25,
    "C": 8.36,
    "H": 6.04,
    "K": 10.45,
    "R": 12.48,
    "Y": 10.46,
}


def charge_at_ph(seq: str, pH: float) -> float:
    charge = 1 / (1 + 10 ** (pH - PKA_VALUES["N_terminus"]))
    charge -= 1 / (1 + 10 ** (PKA_VALUES["C_terminus"] - pH))

    for aa in seq:
        if aa in ["D", "E", "C", "Y"]:
            charge -= 1 / (1 + 10 ** (PKA_VALUES[aa] - pH))
        elif aa in ["K", "R", "H"]:
            charge += 1 / (1 + 10 ** (pH - PKA_VALUES[aa]))

    return charge
