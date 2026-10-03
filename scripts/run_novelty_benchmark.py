import os
import pandas as pd

# The 20 submitted sequences
candidates = [
    'APEX_EGFR_PROTEINBASE_SUBMISSION.csv',
    '/home/agni/ApexEGFR/APEX_EGFR_PROTEINBASE_SUBMISSION.csv',
    os.path.join(os.path.dirname(__file__), '..', 'APEX_EGFR_PROTEINBASE_SUBMISSION.csv')
]
csv_path = next((p for p in candidates if os.path.exists(p)), None)
df_sub = pd.read_csv(csv_path)

# Reference classes
references = {
    'Human EGF': {
        'molecule': 'Endogenous Ligand',
        'pdb': '1IVO, 1JL9',
        'seq': 'NSDSECPLSHDGYCLHDGVCMYIEALDKYACNCVVGYIGERCQYRDLKWWELR'
    },
    'Human TGF-alpha': {
        'molecule': 'Endogenous Ligand',
        'pdb': '1MOX',
        'seq': 'VVSHFNDCPDSHTQFCFHGTCRFLVQEDKPACVCHSGYVGARCEHADLLA'
    },
    'Affibody Z_EGFR:03115': {
        'molecule': 'Engineered 3-Helix Scaffold',
        'pdb': '2KZI',
        'seq': 'VDNKFNKEMWAAWEEIRNLPNLNGWQMTAFIASLVDDPSQSANLLAEAKKLNDAQAPK'
    },
    'Nanobody 7D12': {
        'molecule': 'Single-Domain VHH',
        'pdb': '4KRO, 3G9A',
        'seq': 'QVKLEESGGGSVQTGGSLRLTCAASGRTSRSYGMGWFRQAPGKEREFVSGISWRGDSTGYADSVKGRFTISRDNAKNTVDLQMNSLKPEDTAIYYCAAAAGSAWYGTLYEYDYWGQGTQVTVSS'
    },
    'DARPin E01': {
        'molecule': 'Ankyrin Repeat Protein',
        'pdb': '1SVX, 4AAG',
        'seq': 'DSDLGKKLLEAARAGQDDEVRILMANGADVNATDYWGWTPLHLAAYQGHLEIVEVLLKNGADVNAQDKFGKTAFDISIDNGNEDLAEILQ'
    },
    'Cetuximab CDR-H3': {
        'molecule': 'Therapeutic Monoclonal Ab',
        'pdb': '1YY9',
        'seq': 'ALTYYDYEFAY'
    },
    'Panitumumab CDR-H3': {
        'molecule': 'Therapeutic Monoclonal Ab',
        'pdb': '5SX4',
        'seq': 'DVGYCSSSNCPD'
    },
    'Matuzumab CDR-H3': {
        'molecule': 'Therapeutic Monoclonal Ab',
        'pdb': '3C09',
        'seq': 'QHWYFDL'
    }
}

def max_local_identity(query, target, w=15):
    """Compute maximum percentage sequence identity across sliding windows of length w."""
    win = min(w, len(query), len(target))
    max_id = 0.0
    for i in range(len(query) - win + 1):
        q_sub = query[i:i+win]
        for j in range(len(target) - win + 1):
            t_sub = target[j:j+win]
            matches = sum(1 for a, b in zip(q_sub, t_sub) if a == b)
            ident = matches / win
            if ident > max_id:
                max_id = ident
    return max_id

results = []
for ref_name, ref_data in references.items():
    ref_seq = ref_data['seq']
    max_observed = 0.0
    cdr_matches = 0
    is_cdr = 'CDR-H3' in ref_name
    
    for _, row in df_sub.iterrows():
        s = row['sequence']
        w = len(ref_seq) if is_cdr else 15
        ident = max_local_identity(s, ref_seq, w=w)
        if ident > max_observed:
            max_observed = ident
        if is_cdr and ident >= 0.70:
            cdr_matches += 1
            
    results.append({
        'Reference Target Class': ref_data['molecule'],
        'Representative Molecule': ref_name,
        'PDB Reference': ref_data['pdb'],
        'Max 15-mer Identity': f"{round(max_observed * 100)}%",
        'Match to CDR-H3 Motif': 'Zero match' if is_cdr and cdr_matches == 0 else ('None' if not is_cdr else f'{cdr_matches} matches')
    })

res_df = pd.DataFrame(results)
os.makedirs('outputs', exist_ok=True)
res_df.to_csv('outputs/novelty_benchmark_results.csv', index=False)

print("=== NOVELTY BENCHMARK RESULTS ===")
print(res_df.to_string(index=False))
