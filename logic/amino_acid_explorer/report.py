from data import AMINO_ACID_DATA
from data.amino_acid_explorer.chemical import AMINO_ACID_FORMULAS
from data.amino_acid_explorer.genetics import AA_TO_CODONS_RNA
from data.amino_acid_explorer.structure_images import AMINO_ACID_STRUCTURE_IMAGES

from .blosum62 import get_blosum62_row


def build_amino_acid_report(code: str) -> dict:
    """Build a structured, UI-friendly report for a single amino acid."""
    code = (code or "").upper().strip()
    aa = AMINO_ACID_DATA.get(code)
    if not aa:
        raise KeyError(f"Unknown amino acid code: {code}")

    blosum_row = get_blosum62_row(code)

    return {
        "one_letter": code,
        "three_letter": aa.get("three"),
        "name": aa.get("name"),
        "category": aa.get("category"),
        "molecular_weight": aa.get("mw"),
        "hydropathy_kd": aa.get("kd"),
        "pKa": aa.get("pKa"),
        "formula": AMINO_ACID_FORMULAS.get(code),
        "structure_image": AMINO_ACID_STRUCTURE_IMAGES.get(code),
        "codons_rna": AA_TO_CODONS_RNA.get(code, []),
        "blosum62_row": blosum_row,
    }
