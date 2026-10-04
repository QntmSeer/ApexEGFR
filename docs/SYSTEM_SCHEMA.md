# ApexEGFR System Schema & Data Architecture Specification

This document details the end-to-end computational system architecture, design state machine, and data field schemas for the **ApexEGFR** de novo miniprotein engineering platform.

---

## 1. End-to-End System Pipeline Architecture

```mermaid
flowchart TD
    subgraph S1["1. Target Epitope Preparation"]
        T1["Human EGFR Domain III (PDB: 6ARU)"]
        T2["Mouse EGFR Domain III (AF-Q01279)"]
        H1["Conserved β-Sheet Platform (Leu325, Phe357, Gln384, His409)"]
        T1 --> H1
        T2 --> H1
    end

    subgraph S2["2. Generative Backbone & Sequence Sampling"]
        B1["BindCraft v2 Backbone Scaffolding (68–89 AA)"]
        P1["ProteinMPNN Sequence Sampling (Paratope Optimization)"]
        A1["AlphaFold2 Multimer Co-Folding (Dual-Species + 52-AA Tag)"]
        S1 --> B1 --> P1 --> A1
    end

    subgraph S3["3. Active Learning & Negative Feedback Cloud"]
        F1["Failure Classifier (Clashes, Tag Distortions, Glycan Occlusion)"]
        N1["Negative Density Exclusion Cloud (Paratope Resampling)"]
        A1 --> F1 --> N1
        N1 -- "Feedback Loop" --> P1
    end

    subgraph S4["4. Protonation & pH-Switch Design"]
        PMPNN["Proton-PottsMPNN Microstate Conditioning (HIS-P/HIS-S, ASP-P/ASP-D)"]
        DESEL["ΔE_sel Evaluation (Acidic Tumor pH 6.5 vs Physiological pH 7.4)"]
        A1 --> PMPNN --> DESEL
    end

    subgraph S5["5. All-Atom Molecular Dynamics Telemetry"]
        OMM["OpenMM 8.6 (Amber14SB, TIP3P Water Box, 310 K, 10 ns)"]
        RMSD["RMSD Convergence (1.10 Å) & Contact Persistence (855.8 contacts)"]
        DESEL --> OMM --> RMSD
    end

    subgraph S6["6. Gatekeeper Audit & Sequence Novelty"]
        VERIFY["verify_submission.py (8-Step Gatekeeper Compliance Audit)"]
        TYPER["run_novelty_benchmark.py (ProteinTyper 15-mer Window Audit)"]
        RMSD --> VERIFY
        VERIFY --> TYPER
    end

    subgraph S7["7. Designated Submission & Live Registry"]
        SUB["APEX_EGFR_PROTEINBASE_SUBMISSION.csv (19 Designated Miniproteins)"]
        REG["PROTEINBASE_DESIGNATED_MAPPING.csv (Live Codenames e.g., swift-cat-wave)"]
        REPORT["ApexEGFR_Technical_Report.pdf (29-Page Compiled Whitepaper)"]
        TYPER --> SUB & REG & REPORT
    end

    classDef stage fill:#f9f9f9,stroke:#333,stroke-width:1px;
    classDef highlight fill:#d1fae5,stroke:#059669,stroke-width:2px;
    class S7 highlight;
```

---

## 2. Dataset & CSV Schema Specification

### A. Official Submission Manifest (`APEX_EGFR_PROTEINBASE_SUBMISSION.csv`)

| Field Name | Data Type | Constraint / Format | Description |
| :--- | :--- | :--- | :--- |
| **`name`** | String | `APEX-EGFR-XX` | Unique construct handle submitted to Proteinbase |
| **`sequence`** | String | 68–89 AA (20 Canonical AA) | Full amino acid sequence (Zero free Cysteines, No Tags attached) |

---

### B. Biophysical Audited Metrics Summary (`FINAL_20_METRICS_SUMMARY.csv`)

| Field Name | Data Type | Range / Format | Description |
| :--- | :--- | :--- | :--- |
| **`name`** | String | `APEX-EGFR-XX` | Design identifier |
| **`scaffold_hash`** | String | Hexadecimal (8-char) | Structural backbone fold family identifier |
| **`length`** | Integer | 68–89 AA | Sequence length in amino acids |
| **`iptm_human`** | Float | $0.662\text{--}0.836$ | AF2 predicted interface TM-score against Human EGFR Domain III |
| **`iptm_mouse`** | Float | $0.647\text{--}0.821$ | AF2 predicted interface TM-score against Mouse EGFR Domain III |
| **`tag_clearance_angstrom`** | Float | $4.43\text{--}21.80\text{ \AA}$ | Minimum distance from C-terminal GFP11+TwinStrep tag to target surface |
| **`glycan_clearance_angstrom`** | Float | $15.2\text{--}31.7\text{ \AA}$ | Minimum distance to N-linked glycans ($Asn328$, $Asn420$) |
| **`delta_E_sel`** | Float | $-5.11\text{ to } -6.08\text{ kcal/mol}$ | pH selective binding energy differential (Proton-PottsMPNN) |
| **`novelty_score`** | Float | $0.75\text{--}1.00$ | ProteinTyper 15-mer sliding window novelty score ($\ge 3/4$) |
| **`design_category`** | String | Category Label | Primary biophysical designation (e.g. `Proton-Potts pH-Switch (His-Lead)`) |

---

### C. Live Proteinbase Codename Cross-Reference (`PROTEINBASE_DESIGNATED_MAPPING.csv`)

| Field Name | Data Type | Example Value | Description |
| :--- | :--- | :--- | :--- |
| **`submitted_id`** | String | `APEX-EGFR-01` | Candidate name in competition submission CSV |
| **`proteinbase_codename`** | String | `swift-cat-wave` | Auto-assigned three-word animal/gem/nature codename on Proteinbase |
| **`scaffold_hash`** | String | `4a8595e9` | Fold family scaffold hash |
| **`length`** | Integer | `89` | Amino acid length |
| **`design_category`** | String | `Curated Baseline Lead` | Functional role designation |
