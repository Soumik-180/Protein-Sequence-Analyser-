import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
from pipeline import run_analysis, format_sequence
from core import smith_waterman

st.set_page_config(page_title="Protein Sequence Analyzer", page_icon="🧬", layout="wide")

st.title("Protein Sequence Analyzer")
st.write("A pure Python bioinformatics pipeline. Scientifically accurate algorithms.")

with st.sidebar:
    st.header("Tools")
    mode = st.radio("Select Mode", ["Sequence Analysis", "Sequence Alignment"])

if mode == "Sequence Analysis":
    st.subheader("Input Sequence")
    raw_seq = st.text_area("Enter protein sequence (raw or multi-FASTA)", height=150, 
                           help="Only the first sequence will be analyzed from multi-FASTA inputs.")

    if st.button("Analyze Sequence"):
        if raw_seq:
            with st.spinner("Analyzing..."):
                results = run_analysis(raw_seq)
                if "error" in results:
                    st.error(results["error"])
                else:
                    st.success("Analysis Complete!")
                    
                    st.subheader("Core Biochemical Properties")
                    col1, col2, col3, col4, col5 = st.columns(5)
                    col1.metric("Length", f"{results['length']} aa")
                    col2.metric("Weight", f"{results['molecular_weight']:.2f} Da")
                    col3.metric("pI", f"{results['isoelectric_point']:.2f}")
                    col4.metric("GRAVY", f"{results['gravy']:.2f}")
                    col5.metric("Instability Index", f"{results['instability']['score']:.1f}", 
                                "Stable" if results['instability']['stable'] else "Unstable")
                    
                    c1, c2 = st.columns(2)
                    with c1:
                        st.subheader("Amino Acid Composition")
                        df_comp = pd.DataFrame(list(results['composition'].items()), columns=['AA', 'Percentage (%)'])
                        df_comp = df_comp.sort_values(by='Percentage (%)', ascending=False).reset_index(drop=True)
                        st.dataframe(df_comp, use_container_width=True)
                    
                    with c2:
                        st.subheader("Charge vs pH Curve")
                        curve_data = results['charge_curve']
                        fig_q, ax_q = plt.subplots(figsize=(6, 4))
                        ax_q.plot(curve_data['ph'], curve_data['charge'], color='purple')
                        ax_q.axhline(0, color='gray', linestyle='--')
                        ax_q.axvline(results['isoelectric_point'], color='red', linestyle='--', label=f"pI = {results['isoelectric_point']:.2f}")
                        ax_q.set_xlabel("pH")
                        ax_q.set_ylabel("Net Charge")
                        ax_q.legend()
                        st.pyplot(fig_q)

                    st.subheader("Hydrophobicity (Kyte-Doolittle)")
                    if results['hydrophobicity_plot']:
                        df_hydro = pd.DataFrame(results['hydrophobicity_plot'])
                        fig, ax = plt.subplots(figsize=(10, 3))
                        ax.fill_between(df_hydro['position'], df_hydro['score'], 0, where=(df_hydro['score'] > 0), color='salmon', alpha=0.5, interpolate=True)
                        ax.fill_between(df_hydro['position'], df_hydro['score'], 0, where=(df_hydro['score'] <= 0), color='lightblue', alpha=0.5, interpolate=True)
                        ax.plot(df_hydro['position'], df_hydro['score'], color='black', linewidth=1)
                        ax.axhline(0, color='gray', linestyle='--')
                        ax.set_xlabel("Position")
                        ax.set_ylabel("Hydropathy Score")
                        st.pyplot(fig)
                        
                    c3, c4 = st.columns(2)
                    with c3:
                        st.subheader("Detected Motifs")
                        if results['motifs']:
                            for m in results['motifs']:
                                st.write(f"- **{m['type']}** at pos {m['position']}")
                        else:
                            st.write("No common motifs detected.")
                    with c4:
                        st.subheader("Secondary Structure (Proxy)")
                        if results['secondary_structure']:
                            structs = results['secondary_structure']
                            alpha = structs.count('Alpha Helix') / len(structs) * 100
                            beta = structs.count('Beta Sheet') / len(structs) * 100
                            coil = structs.count('Coil') / len(structs) * 100
                            st.write(f"- **Alpha Helix:** {alpha:.1f}%")
                            st.write(f"- **Beta Sheet:** {beta:.1f}%")
                            st.write(f"- **Coil:** {coil:.1f}%")

elif mode == "Sequence Alignment":
    st.subheader("🧬 Smith-Waterman Local Alignment")
    seq1 = st.text_area("Sequence 1 (raw or FASTA)")
    seq2 = st.text_area("Sequence 2 (raw or FASTA)")
    
    if st.button("Align Sequences"):
        valid_aa = set("ACDEFGHIKLMNPQRSTVWY")
        s1 = format_sequence(seq1)
        s2 = format_sequence(seq2)
        
        if not s1 or not s2:
            st.warning("Please provide valid sequences in both fields.")
        elif not all(aa in valid_aa for aa in s1) or not all(aa in valid_aa for aa in s2):
            st.error("Invalid amino acids detected in your sequence.")
        else:
            alignment = smith_waterman(s1, s2)
            st.success("Alignment Complete!")
            col1, col2, col3 = st.columns(3)
            col1.metric("Alignment Score", alignment['score'])
            col2.metric("Similarity", f"{alignment['similarity']:.1f}%")
            col3.metric("Matches", f"{alignment['matches']} / {alignment['alignLen']}")

            st.write("### Alignment Trace")
            if alignment['alignLen'] > 0:
                match_str = ""
                for a, b in zip(alignment['aligned_seq1'], alignment['aligned_seq2']):
                    if a == b: match_str += "|"
                    elif a == '-' or b == '-': match_str += " "
                    else: match_str += "."
                st.text(f"Seq1: {alignment['aligned_seq1']}\n      {match_str}\nSeq2: {alignment['aligned_seq2']}")
            else:
                st.write("No positive-scoring alignment path found.")
