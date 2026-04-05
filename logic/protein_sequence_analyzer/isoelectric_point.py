from .charge_model import charge_at_ph


def calculate_pi(seq: str) -> float:
    low_ph, high_ph, pi = 0.0, 14.0, 7.0
    for _ in range(100):
        pi = (low_ph + high_ph) / 2
        if charge_at_ph(seq, pi) > 0:
            low_ph = pi
        else:
            high_ph = pi
    return pi
