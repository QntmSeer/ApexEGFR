import os
import pandas as pd

CODENAME_MAP = {
    'APEX-EGFR-01': 'swift-cat-wave',
    'APEX-EGFR-02': 'soft-goat-topaz',
    'APEX-EGFR-03': 'strong-ram-opal',
    'APEX-EGFR-04': 'silver-eagle-ruby',
    'APEX-EGFR-05': 'steady-heron-opal',
    'APEX-EGFR-06': 'dark-ant-reed',
    'APEX-EGFR-08': 'wild-owl-fern',
    'APEX-EGFR-09': 'rough-tiger-ruby',
    'APEX-EGFR-10': 'crimson-owl-orchid',
    'APEX-EGFR-11': 'rough-ox-birch',
    'APEX-EGFR-12': 'wild-otter-maple',
    'APEX-EGFR-13': 'pale-goat-moss',
    'APEX-EGFR-14': 'gentle-crow-dust',
    'APEX-EGFR-15': 'crimson-gecko-birch',
    'APEX-EGFR-16': 'quick-deer-ruby',
    'APEX-EGFR-17': 'scarlet-ant-ash',
    'APEX-EGFR-18': 'crimson-ram-granite',
    'APEX-EGFR-19': 'green-ram-oak',
    'APEX-EGFR-20': 'vast-eagle-lotus'
}

def main():
    summary_path = os.path.join('data', 'FINAL_20_METRICS_SUMMARY.csv')
    if not os.path.exists(summary_path):
        summary_path = 'FINAL_20_METRICS_SUMMARY.csv'
        
    df = pd.read_csv(summary_path)
    df_19 = df[df['submission_name'] != 'APEX-EGFR-07'].copy().reset_index(drop=True)
    df_19['proteinbase_codename'] = df_19['submission_name'].map(CODENAME_MAP)
    
    # Save CSV
    cols = ['submission_name', 'proteinbase_codename', 'scaffold', 'length', 'category',
            'human_tagged_iptm', 'mouse_tagged_iptm', 'tag_to_epitope_dist_A', 'sequence']
    df_out = df_19[cols].copy()
    df_out.columns = ['submission_id', 'proteinbase_codename', 'scaffold', 'length', 'category',
                      'human_tagged_iptm', 'mouse_tagged_iptm', 'tag_dist_A', 'sequence']
    
    out_csv = os.path.join('data', 'PROTEINBASE_DESIGNATED_MAPPING.csv')
    df_out.to_csv(out_csv, index=False)
    print(f"Wrote {out_csv}")

    # Generate Markdown documentation
    md = [
        "# Official Proteinbase Registry & Codename Cross-Reference",
        "",
        "**Author**: `QntmSeer`  ",
        "**Collection**: `[Anthropic × Adaptyv Protein Design Competition] Submission 1`  ",
        "**Designated Status**: **DESIGNATED (19 Verified Constructs)**  ",
        "",
        "This registry maps our internal campaign identifiers (`APEX-EGFR-XX`) to Adaptyv Bio's auto-assigned Proteinbase slugs (`adjective-animal-gem/mineral/plant`) visible on the live platform.",
        "",
        "| Submitted ID | Live Proteinbase Codename | Scaffold | Length | Design Class | Human $i\\_pTM$ | Mouse $i\\_pTM$ | Tag Clearance |",
        "| :--- | :--- | :---: | :---: | :--- | :---: | :---: | :---: |"
    ]
    
    for _, r in df_out.iterrows():
        md.append(f"| **`{r['submission_id']}`** | `{r['proteinbase_codename']}` | `{r['scaffold']}` | {r['length']} AA | {r['category']} | {r['human_tagged_iptm']:.3f} | {r['mouse_tagged_iptm']:.3f} | {r['tag_dist_A']:.1f} Å |")
    
    md.append("")
    md.append("### Key Lead Cross-Reference")
    md.append("- **10 ns Explicit-Solvent MD Lead**: `APEX-EGFR-16` $\\leftrightarrow$ **`quick-deer-ruby`** (Scaffold `950df084`, 68 AA, $i\\_pTM = 0.836 / 0.821$)")
    md.append("- **Proton-Potts His-Switch Lead**: `APEX-EGFR-02` $\\leftrightarrow$ **`soft-goat-topaz`** (Scaffold `89062c71`, 79 AA, $\\Delta E_{\\text{sel}} = -5.11$)")
    md.append("- **Proton-Potts Asp-Switch Lead**: `APEX-EGFR-06` $\\leftrightarrow$ **`dark-ant-reed`** (Scaffold `89062c71`, 79 AA, $\\Delta E_{\\text{sel}} = -6.08$)")
    md.append("- **Ultra-High Tag Clearance Controls**: `APEX-EGFR-19` $\\leftrightarrow$ **`green-ram-oak`** & `APEX-EGFR-20` $\\leftrightarrow$ **`vast-eagle-lotus`** (Scaffold `991afb06`, 89 AA, $21.8\\text{ Å}$ tag buffer, 0 histidines)")
    md.append("")
    
    out_md = os.path.join('docs', 'PROTEINBASE_DESIGNATED_MAPPING.md')
    with open(out_md, 'w', encoding='utf-8') as f:
        f.write("\n".join(md))
    print(f"Wrote {out_md}")

if __name__ == '__main__':
    main()
