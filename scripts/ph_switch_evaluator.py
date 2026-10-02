"""
ph_switch_evaluator.py — Biophysical Scoring of pH 6.5 vs 7.4 Binding Switch
=============================================================================
Calculates the predicted binding free energy difference (Delta-Delta-G_pH)
between acidic tumor microenvironment (pH 6.5) and neutral healthy tissue (pH 7.4).

Key biophysical determinants:
1. Interfacial Histidines within 4.5 A of target acidic residues (Glu472, Asp355, Glu370).
2. Shifted pKa of interface Histidines when in proximity to negative carboxylates (pKa ~ 6.6 - 6.8).
3. Fractional protonation:
   - At pH 6.5: ~60-75% protonated (His+) -> active salt bridge formation
   - At pH 7.4: ~10-15% protonated (His0) -> loss of salt bridge, desolvation penalty
4. Predicted Delta-Delta-G_pH = Delta-G_bind(pH 6.5) - Delta-G_bind(pH 7.4)
   Target: Delta-Delta-G_pH <= -2.0 kcal/mol (indicating >25-fold affinity enhancement at pH 6.5)
"""

import math
import numpy as np

def evaluate_ph_switch(binder_coords, binder_resnames, target_coords, target_resnames, target_resnums):
    """
    Evaluates interface contacts in a complex structure:
    binder_coords: (N, 3) numpy array
    binder_resnames: list of 3-letter codes for binder residues
    target_coords: (M, 3) numpy array
    target_resnames: list of 3-letter codes for target residues
    target_resnums: list of 1-based target residue numbers
    """
    # Identify target acidic anchors (specifically Glu472 on Domain III, plus Asp/Glu in contact zone)
    acidic_contacts = []
    
    # Distance matrix between binder atoms and target atoms
    for i, (b_coord, b_res) in enumerate(zip(binder_coords, binder_resnames)):
        if b_res != 'HIS':
            continue
        # Find closest target atom
        dists = np.linalg.norm(target_coords - b_coord, axis=1)
        min_idx = np.argmin(dists)
        min_d = dists[min_idx]
        
        t_res = target_resnames[min_idx]
        t_num = target_resnums[min_idx]
        
        if min_d <= 5.0 and t_res in ['GLU', 'ASP']:
            acidic_contacts.append({
                "binder_atom_idx": i,
                "target_res": f"{t_res}{t_num}",
                "distance": min_d,
                "is_glu472": (t_num == 472)
            })
            
    # Calculate score
    # Each His-Asp/Glu pair at interface contributes:
    # At pH 6.5: favorable salt bridge ~ -3.5 kcal/mol (scaled by distance)
    # At pH 7.4: neutral His with no salt bridge ~ -0.5 kcal/mol
    ddg_total = 0.0
    for ac in acidic_contacts:
        d = ac["distance"]
        dist_factor = max(0.2, min(1.0, (5.0 - d) / 2.0))
        # Salt bridge bonus at pH 6.5 vs pH 7.4
        weight = 1.5 if ac["is_glu472"] else 1.0
        ddg = -3.0 * dist_factor * weight
        ddg_total += ddg
        
    num_his_acidic_pairs = len(acidic_contacts)
    is_switchable = (ddg_total <= -2.0) and (num_his_acidic_pairs >= 1)
    
    return {
        "is_switchable": is_switchable,
        "predicted_ddg_ph": round(ddg_total, 2),
        "num_his_acidic_pairs": num_his_acidic_pairs,
        "details": acidic_contacts
    }

if __name__ == "__main__":
    print("pH Switch Evaluator module ready.")
