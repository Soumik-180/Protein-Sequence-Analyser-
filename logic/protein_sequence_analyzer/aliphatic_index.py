def calculate_aliphatic_index(seq: str) -> float:
    if not seq:
        return 0

    valid_seq = [aa for aa in seq if aa in set("ACDEFGHIKLMNPQRSTVWY")]
    if not valid_seq:
        return 0

    a_count = valid_seq.count("A")
    v_count = valid_seq.count("V")
    i_count = valid_seq.count("I")
    l_count = valid_seq.count("L")

    return ((a_count + 2.9 * v_count + 3.9 * (i_count + l_count)) / len(valid_seq)) * 100
