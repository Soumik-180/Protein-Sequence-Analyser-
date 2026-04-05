from .charge_model import charge_at_ph


def compute_charge_curve(seq: str) -> dict[str, list[float]]:
    ph_values = [i * 0.5 for i in range(0, 29)]
    charges = [charge_at_ph(seq, ph) for ph in ph_values]
    return {"ph": ph_values, "charge": charges}
