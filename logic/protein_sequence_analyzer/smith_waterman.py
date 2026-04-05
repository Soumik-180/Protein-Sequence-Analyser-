from data import BLOSUM62


def smith_waterman(seq1: str, seq2: str, gap_open: int = -11, gap_ext: int = -1) -> dict:
    n, m = len(seq1), len(seq2)
    if n == 0 or m == 0:
        return {
            "score": 0,
            "similarity": 0,
            "matches": 0,
            "alignLen": 0,
            "aligned_seq1": "",
            "aligned_seq2": "",
        }

    H = [[0] * (m + 1) for _ in range(n + 1)]
    E = [[0] * (m + 1) for _ in range(n + 1)]
    F = [[0] * (m + 1) for _ in range(n + 1)]

    # State flags: 0=Stop, 1=H(diag), 2=E(up), 3=F(left)
    TB_H = [[0] * (m + 1) for _ in range(n + 1)]
    TB_E = [[0] * (m + 1) for _ in range(n + 1)]
    TB_F = [[0] * (m + 1) for _ in range(n + 1)]

    max_score = 0
    max_i, max_j = 0, 0

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            match_score = BLOSUM62.get(seq1[i - 1], {}).get(seq2[j - 1], -4)

            # E (gap in seq2 / vertical)
            e_open = H[i - 1][j] + gap_open
            e_ext = E[i - 1][j] + gap_ext
            if e_open >= e_ext:
                E[i][j] = e_open
                TB_E[i][j] = 1  # from H
            else:
                E[i][j] = e_ext
                TB_E[i][j] = 2  # from E

            # F (gap in seq1 / horizontal)
            f_open = H[i][j - 1] + gap_open
            f_ext = F[i][j - 1] + gap_ext
            if f_open >= f_ext:
                F[i][j] = f_open
                TB_F[i][j] = 1  # from H
            else:
                F[i][j] = f_ext
                TB_F[i][j] = 3  # from F

            # H
            diag = H[i - 1][j - 1] + match_score
            scores = (0, diag, E[i][j], F[i][j])
            best = max(scores)
            H[i][j] = best

            if best == 0:
                TB_H[i][j] = 0
            elif best == diag:
                TB_H[i][j] = 1
            elif best == E[i][j]:
                TB_H[i][j] = 2
            elif best == F[i][j]:
                TB_H[i][j] = 3

            if H[i][j] > max_score:
                max_score = H[i][j]
                max_i, max_j = i, j

    i, j = max_i, max_j
    matches, align_len = 0, 0
    al1, al2 = [], []
    state = 1  # Start traceback in H state

    while i > 0 and j > 0:
        if state == 1:  # In H matrix
            if TB_H[i][j] == 0:
                break
            tr = TB_H[i][j]
            if tr == 1:
                align_len += 1
                if seq1[i - 1] == seq2[j - 1]:
                    matches += 1
                al1.append(seq1[i - 1])
                al2.append(seq2[j - 1])
                i -= 1
                j -= 1
            elif tr == 2:
                state = 2
            elif tr == 3:
                state = 3
        elif state == 2:  # In E matrix
            align_len += 1
            al1.append(seq1[i - 1])
            al2.append("-")
            state = TB_E[i][j]
            i -= 1
        elif state == 3:  # In F matrix
            align_len += 1
            al1.append("-")
            al2.append(seq2[j - 1])
            state = TB_F[i][j]
            j -= 1

    similarity = (matches / align_len * 100) if align_len > 0 else 0
    return {
        "score": max_score,
        "similarity": similarity,
        "matches": matches,
        "alignLen": align_len,
        "aligned_seq1": "".join(reversed(al1)),
        "aligned_seq2": "".join(reversed(al2)),
    }
