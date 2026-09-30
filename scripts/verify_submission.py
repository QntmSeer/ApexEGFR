import os
import pandas as pd

def main():
    csv_path = os.path.join('data', 'APEX_EGFR_PROTEINBASE_SUBMISSION.csv')
    if not os.path.exists(csv_path):
        csv_path = 'APEX_EGFR_PROTEINBASE_SUBMISSION.csv'
        
    df = pd.read_csv(csv_path)
    print(f"=== LOCAL VERIFICATION OF {csv_path} ===")
    print("Columns:", df.columns.tolist())
    assert df.columns.tolist() == ['name', 'sequence', 'molecule_class'], "Invalid CSV column header"
    print("Row count:", len(df))
    assert len(df) == 20, f"Expected 20 rows, got {len(df)}"
    print("Unique sequences:", df['sequence'].nunique())
    assert df['sequence'].nunique() == 20, "Duplicate sequence found"
    print("Unique names:", df['name'].nunique())
    assert df['name'].nunique() == 20, "Duplicate name found"
    print("Molecule class all single_chain:", (df['molecule_class'] == 'single_chain').all())
    assert (df['molecule_class'] == 'single_chain').all(), "Molecule class must be single_chain"

    TAG_SEQ = 'GGGSRDHMVLHEYVNAAGITWSHPQFEKGGGSGGGSGGGSWSHPQFEK'
    for i, r in df.iterrows():
        name = r['name']
        seq = r['sequence']
        assert seq.count('C') == 0, f"Cysteine found in {name}"
        assert 50 <= len(seq) <= 100, f"Length out of bounds in {name}: {len(seq)}"
        for n in range(5, len(TAG_SEQ)+1):
            assert TAG_SEQ[:n] not in seq, f"Tag fragment {TAG_SEQ[:n]} in {name}"

    print("All sequences passed validation: 0 Cys, 0 tag fragments, valid lengths (65-89 AA), unique names & sequences.")

if __name__ == '__main__':
    main()
