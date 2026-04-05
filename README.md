# Protein Sequence Analyzer

A modular, pure-Python bioinformatics web application built with Streamlit. It performs rigorous physicochemical analysis, sequence characterization, and local alignment for protein sequences with scientifically grounded methods.

---

## 🔗 Live Demo

Access the deployed application here:  
https://protein-sequence-analyser.streamlit.app/

---

## Features

### 1. Physicochemical Analysis
- **Molecular Weight:** Correctly computed with peptide bond correction (water loss).
- **Isoelectric Point (pI):** Estimated using Henderson–Hasselbalch equilibrium with iterative charge balancing.
- **Charge Curve:** Visualizes net charge across pH 0–14.
- **GRAVY Score:** Based on the Kyte–Doolittle hydropathy scale.
- **Instability Index:** Calculated using the Guruprasad et al. (1990) DIWV dipeptide matrix.
- **Aliphatic Index:** Reflects thermostability based on aliphatic residue composition.

---

### 2. Sequence Analysis & Visualization
- **Hydrophobicity Plot:** Sliding-window Kyte–Doolittle profile (default window = 9).
- **Amino Acid Composition:** Residue counts and percentage distribution.
- **Motif Detection:** Regex-based detection of common motifs (e.g., N-glycosylation, phosphorylation sites).
- **Secondary Structure (Heuristic):**
  - Rule-based approximation using residue propensities
  - **Note:** This is NOT a full Chou–Fasman or ML-based predictor.

---

### 3. Local Sequence Alignment
- **Smith–Waterman Algorithm:**
  - Local alignment with affine gap penalties
  - Uses three matrices: H (score), E (gap in seq1), F (gap in seq2)
- **BLOSUM62 Matrix:** Standard substitution scoring
- **Alignment Output:**
  - Aligned sequences
  - Match/mismatch visualization
  - Alignment score and identity %

---

### 4. Amino Acid Explorer
- **Single Amino Acid Lookup:** Accepts 1-letter, 3-letter, or full name
- **Chemical Reference:** Formula (with subscripts) + 2D structure image (from `structure/`)
- **Genetics:** RNA codons table
- **BLOSUM62 Explorer:** Interactive substitution score comparison

---

## Architecture

The application follows a clean, modular pipeline design:

```
Input → Parsing/Validation (pipeline) → Computation (logic + data) → Results → UI Rendering
```

### Components:
- **`app.py`** → Streamlit UI layer
- **`pipeline.py`** → Input validation, FASTA parsing, orchestration
- **`logic/`** → Scientific computation modules (alignment, pI, hydrophobicity, etc.)
- **`data/`** → Static biological datasets (AA properties, BLOSUM62, DIWV)

### Module Map

**`logic/protein_sequence_analyzer/`** (one file per calculation)
- `molecular_weight.py` → `calculate_molecular_weight`
- `composition.py` → `calculate_composition`
- `charge_model.py` → `charge_at_ph` (shared)
- `charge_curve.py` → `compute_charge_curve`
- `isoelectric_point.py` → `calculate_pi`
- `gravy.py` → `calculate_gravy`
- `instability.py` → `calculate_instability_index`
- `aliphatic_index.py` → `calculate_aliphatic_index`
- `hydrophobicity.py` → `calculate_hydrophobicity`
- `motifs.py` → `detect_motifs`
- `secondary_structure.py` → `simple_secondary_structure_prediction`
- `smith_waterman.py` → `smith_waterman`

**`logic/amino_acid_explorer/`** (one file per helper)
- `normalize.py` → `normalize_amino_acid_query`
- `report.py` → `build_amino_acid_report`
- `blosum62.py` → `get_blosum62_row`, `get_blosum62_score`

**`data/protein_sequence_analyzer/`** (datasets)
- `amino_acids.py` → `AMINO_ACID_DATA`
- `blosum62.py` → `BLOSUM62`
- `diwv.py` → `DIWV`

**`data/amino_acid_explorer/`** (offline reference tables)
- `chemical.py` → `AMINO_ACID_FORMULAS`
- `structure_images.py` → `AMINO_ACID_STRUCTURE_IMAGES`
- `genetics.py` → `AA_TO_CODONS_RNA`

**`data/`**
- `__init__.py` re-exports `AMINO_ACID_DATA`, `BLOSUM62`, `DIWV`

---

## Input Handling

- Accepts raw sequences or FASTA format
- Automatically strips headers (`>...`) and whitespace
- Validates amino acid characters (A–Z, standard residues only)
- Returns structured errors for invalid input

---

## 🚀 Setup

### Requirements
- Python 3.9+

### Installation

```bash
git clone <your-repo-url>
cd Protein
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

---

## Run the App

```bash
streamlit run app.py
```

Open: http://localhost:8501

---

## Output

The tool provides:
- Numerical metrics (MW, pI, GRAVY, instability, aliphatic index)
- Graphs (hydrophobicity, charge vs pH)
- Alignment results (score + formatted alignment)
- Motif annotations

---

## Limitations

- Secondary structure prediction is heuristic (not Chou–Fasman or ML-based)
- No 3D structure modeling
- No database integration (e.g., UniProt, PDB)

---

## Future Improvements

- Full Chou–Fasman or ML-based secondary structure prediction
- FASTA multi-sequence support
- Export (CSV / JSON reports)
- REST API backend (FastAPI)
- Integration with external databases (UniProt, PDB)

---

## License

For academic and educational use.

Developed by **Soumik Ray**  
B.Tech Biotechnology (Computational Biology)  
Sharda University