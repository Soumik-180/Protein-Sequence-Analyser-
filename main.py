from pipeline import run_analysis

if __name__ == "__main__":
    test_seq = """>sp|P01112|HRAS_HUMAN GTPase HRas OS=Homo sapiens
    MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAG
    QEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHQYREQIKRVKDSDDVPMVLVGNKCDL
    AARTVESRQAQDLARSYGIPYIETSAKTRQGVEDAFYTLVREIRQHKLRKLNPPDESGPG
    CMSCKCVLS
    """
    results = run_analysis(test_seq)
    
    print("--- Basic Stats ---")
    print(f"Length: {results.get('length', 'Error')}")
    print(f"Molecular Weight: {results.get('molecular_weight', 'Error')}")
    print(f"Isoelectric Point: {results.get('isoelectric_point', 'Error')}")
    
    print("\n--- Indices ---")
    print(f"GRAVY: {results.get('gravy', 'Error')}")
    print(f"Aliphatic Index: {results.get('aliphatic_index', 'Error')}")
    
    instability = results.get('instability')
    if instability:
        print(f"Instability Index: {instability['score']:.2f}")
        print(f"Stable: {instability['stable']}")
