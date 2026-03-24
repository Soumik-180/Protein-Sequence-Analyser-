# Protein Sequence Analyzer

A comprehensive, pure-Python bioinformatics web application built with Streamlit. This tool provides deep computational analysis, physicochemical property calculation, and local sequence alignment for protein sequences.

## 🧬 Features

### 1. Physicochemical Analysis
- **Molecular Weight Calculation:** Accurately computes the molecular weight of the protein sequence, adjusting for water loss during peptide bond formation.
- **Isoelectric Point (pI) & Charge Curve:** Calculates the pI and plots the net charge of the protein across the standard pH range (0 to 14) using empirical pKa values.
- **GRAVY Score:** Calculates the Grand Average of Hydropathy using the strict Kyte-Doolittle (`kd`) scale.
- **Instability Index:** Evaluates sequence stability based on the Guruprasad et al. (1990) DIWV (Dipeptide Instability Weight Value) 400-item pairwise matrix.
- **Aliphatic Index:** Computes the relative volume occupied by aliphatic side chains (Alanine, Valine, Isoleucine, and Leucine).

### 2. Sequence Visualization & Motifs
- **Hydrophobicity Plot:** Generates a sliding-window (default window size 9) hydrophobicity chart along the sequence length.
- **Motif Detection:** Scans the sequence for conserved biological motifs using regular expressions (e.g., N-glycosylation sites, Casein kinase II phosphorylation sites).
- **Secondary Structure Prediction:** Provides a heuristic prediction for local secondary structures (Alpha Helix, Beta Sheet, or Coil) across the sequence.

### 3. Local Sequence Alignment
- **Smith-Waterman Algorithm:** Performs local sequence alignment with strict **Affine Gap Penalties**, using 3 independent dynamic programming tracking matrices (H, E, and F).
- **BLOSUM62 Matrix:** Evaluates amino acid substitutions utilizing the standard BLOSUM62 scoring matrix.
- **Sequence Pre-processing:** Automatically strips FASTA headers (`>...`) and whitespace to deliver clean alignments.

## 🏗️ Project Architecture

The codebase is logically partitioned to separate UI, data, and mathematical logic:

* **`app.py`**: The Streamlit user interface. Responsible for rendering inputs, managing state, and visualizing plots/data charts.
* **`pipeline.py`**: The middle tier. Handles validation, FASTA string stripping, error processing, and routes commands between the UI and computation core.
* **`core.py`**: The computation engine. Contains the raw mathematical algorithms (Smith-Waterman backtracking, pH iterators, sequence window aggregation).
* **`data.py`**: The biological data repository. Houses structural definitions like `AMINO_ACID_DATA` (weights, pKa, Kyte-Doolittle scaling), the `BLOSUM62` matrix, and the massive 20x20 `DIWV` instability index matrix.

## 🚀 Setup and Installation

### Prerequisites
- Python 3.9+ (Works flawlessly on macOS / Python 3.14 via lazy bounds mapping).

### Installation Instructions

1. **Clone the repository / navigate to the folder:**
   ```bash
   cd /Users/soumikray/Documents/Protein
   ```

2. **Create and activate a virtual environment:**
   ```bash
   python3 -m venv venv
   source venv/bin/activate
   ```

3. **Install the dependencies:**
   ```bash
   pip install -r requirements.txt
   ```
   *(Note: The requirements utilize `>=` bounds to accommodate modern Python environments and ensure `pandas`, `numpy`, `matplotlib`, and `streamlit` build safely).*

## 💻 Running the Application

To launch the web interface locally, ensure your virtual environment is activated and run:

```bash
streamlit run app.py
```

Streamlit will boot up a local server and give you a `localhost` URL (usually `http://localhost:8501`) that automatically opens in your web browser.

## 📝 Recent Technical Milestones

- **Strict Scaling:** Removed arbitrary biological scalar metrics in favor of pure scientific standards (e.g., exclusively utilizing the Kyte-Doolittle hydropathy scale).
- **Corrected Traceback Matrix:** Advanced DP states within Smith-Waterman natively trace accurate `align_len` offsets without double-counting gap extensions.
- **Data Encapsulation:** Decoupled `DIWV` matrices entirely from the execution codebase to keep operational files lightweight and clean.

# Protein Sequence Analyzer

A modular, pure-Python bioinformatics web application built with Streamlit. This tool performs rigorous computational analysis, physicochemical characterization, and local sequence alignment for protein sequences.

---

## 🧬 Features

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

## 🏗️ Architecture

The application follows a clean, modular pipeline design:

```
Input → Validation → Core Computation → Pipeline → UI Rendering
```

### Components:
- **`app.py`** → Streamlit UI layer
- **`pipeline.py`** → Input validation, FASTA parsing, orchestration
- **`core.py`** → Scientific computation (alignment, pI, hydrophobicity, etc.)
- **`data.py`** → Static biological datasets (AA properties, BLOSUM62, DIWV)

---

## ⚙️ Input Handling

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

## ▶️ Run the App

```bash
streamlit run app.py
```

Open: http://localhost:8501

---

## 📊 Output

The tool provides:
- Numerical metrics (MW, pI, GRAVY, instability, aliphatic index)
- Graphs (hydrophobicity, charge vs pH)
- Alignment results (score + formatted alignment)
- Motif annotations

---

## ⚠️ Limitations

- Secondary structure prediction is heuristic (not Chou–Fasman or ML-based)
- No 3D structure modeling
- No database integration (e.g., UniProt, PDB)

---

## 🧠 Future Improvements

- Full Chou–Fasman or ML-based secondary structure prediction
- FASTA multi-sequence support
- Export (CSV / JSON reports)
- REST API backend (FastAPI)
- Integration with external databases (UniProt, PDB)

---

## 🧾 License

For academic and educational use.