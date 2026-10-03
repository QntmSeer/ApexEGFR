import pandas as pd
import numpy as np
import subprocess
import os

def generate_plots():
    print("=== Generating Data for Gnuplot ===")
    os.makedirs("plots", exist_ok=True)
    df = pd.read_csv("results_backup/all_refolded_consolidated.csv")
    
    def parse_duo_float(val, idx=0):
        try:
            return float(str(val).split(';')[idx])
        except:
            return np.nan

    df['iptm_human'] = df['i_pTM'].apply(lambda x: parse_duo_float(x, 0))
    df['iptm_mouse'] = df['i_pTM'].apply(lambda x: parse_duo_float(x, 1))
    df['buried_area_human'] = df['Interface_BuriedArea'].apply(lambda x: parse_duo_float(x, 0))
    df['is_pass'] = (df['outcome'] == 'passed').astype(int)
    
    # 1. Cross-species data file
    with open("plots/cross_species_pass.dat", "w") as fp, open("plots/cross_species_fail.dat", "w") as ff:
        fp.write("# iptm_human iptm_mouse buried_area\n")
        ff.write("# iptm_human iptm_mouse buried_area\n")
        for _, r in df.dropna(subset=['iptm_human', 'iptm_mouse']).iterrows():
            line = f"{r['iptm_human']:.3f} {r['iptm_mouse']:.3f} {r['buried_area_human']:.1f}\n"
            if r['is_pass'] == 1:
                fp.write(line)
            else:
                ff.write(line)
                
    # 2. AA Enrichment data file
    from collections import Counter
    passed = df[df['outcome'] == 'passed']
    rejected = df[df['outcome'] == 'rejected']
    passed_counts = Counter("".join(passed['Binder_Sequence'].dropna()))
    rejected_counts = Counter("".join(rejected['Binder_Sequence'].dropna()))
    total_passed = sum(passed_counts.values())
    total_rejected = sum(rejected_counts.values())
    
    standard_aas = list("ACDEFGHIKLMNPQRSTVWY")
    aa_data = []
    for aa in standard_aas:
        fp_f = passed_counts[aa] / total_passed
        ff_f = rejected_counts[aa] / total_rejected
        lor = np.log2((fp_f + 1e-4) / (ff_f + 1e-4))
        aa_data.append((aa, lor))
    aa_data.sort(key=lambda x: x[1]) # ascending
    
    with open("plots/aa_enrichment.dat", "w") as f:
        f.write("# index AA log2_odds\n")
        for idx, (aa, lor) in enumerate(aa_data):
            f.write(f"{idx} {aa} {lor:.3f}\n")
            
    # 3. Buried SASA vs Human i_pTM
    with open("plots/sasa_iptm_pass.dat", "w") as fp, open("plots/sasa_iptm_fail.dat", "w") as ff:
        fp.write("# buried_area iptm_human\n")
        ff.write("# buried_area iptm_human\n")
        for _, r in df.dropna(subset=['buried_area_human', 'iptm_human']).iterrows():
            line = f"{r['buried_area_human']:.1f} {r['iptm_human']:.3f}\n"
            if r['is_pass'] == 1:
                fp.write(line)
            else:
                ff.write(line)
                
    print("Local dat files created in plots/")

if __name__ == "__main__":
    generate_plots()
