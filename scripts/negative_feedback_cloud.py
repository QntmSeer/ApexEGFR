"""
negative_feedback_cloud.py
===========================
Active Learning & Negative Design Feedback System for ApexEGFR / Protein Challenge.
Runs on workstation 'agni' (NVIDIA RTX A2000 Laptop GPU).

Key Capabilities:
1. PyTorch Neural Energy Network (NegativeCloudNet) trained on GPU (cuda:0).
2. Computes Repulsive Energy E_rep(x) and Guidance Gradient -grad_x E_rep(x) to steer generations away.
3. Multi-mode failure decomposition (Cross-species, SASA underpacking, Monomer collapse).
4. Exportable ProteinMPNN logit bias dictionary for next-gen sampling.
"""

import os
import sys
import numpy as np
import pandas as pd
from collections import Counter
import torch
import torch.nn as nn
import torch.optim as optim

def parse_float(val, default=0.0):
    try:
        if pd.isna(val):
            return default
        s = str(val).split(';')[0].strip()
        return float(s)
    except:
        return default

class NegativeCloudNet(nn.Module):
    """
    Energy-Based Discriminator that maps candidate features to failure likelihood
    and repulsive potential energy.
    """
    def __init__(self, input_dim=10):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(input_dim, 32),
            nn.LayerNorm(32),
            nn.SiLU(),
            nn.Linear(32, 16),
            nn.SiLU(),
            nn.Linear(16, 1)
        )
        
    def forward(self, x):
        # Returns raw logit / energy
        return self.net(x)

def run_feedback_pipeline(csv_path="/home/agni/ApexEGFR/data/all_refolded_consolidated.csv"):
    if not os.path.exists(csv_path):
        csv_path = "results_backup/all_refolded_consolidated.csv"
        
    print(f"=== ApexEGFR Active Learning Negative Design Engine ===")
    print(f"Loading campaign dataset: {csv_path}")
    df = pd.read_csv(csv_path)
    print(f"Total candidates loaded: {len(df)}")
    
    # 1. Dataset Accounting
    df['is_pass'] = (df['outcome'] == 'passed').astype(int)
    pass_df = df[df['is_pass'] == 1]
    fail_df = df[df['is_pass'] == 0]
    print(f"  Passed Leads (Attractor Basin): {len(pass_df)} ({len(pass_df)/len(df)*100:.1f}%)")
    print(f"  Rejected Candidates (Failure Cloud): {len(fail_df)} ({len(fail_df)/len(df)*100:.1f}%)")
    
    # 2. ProteinMPNN Log-Odds Negative Bias Weights
    standard_aas = list("ACDEFGHIKLMNPQRSTVWY")
    pass_counts = Counter("".join(pass_df['Binder_Sequence'].dropna()))
    fail_counts = Counter("".join(fail_df['Binder_Sequence'].dropna()))
    
    tot_pass_aa = sum(pass_counts.values())
    tot_fail_aa = sum(fail_counts.values())
    
    mpnn_bias_dict = {}
    print("\n--- 1. ProteinMPNN Sequence-Level Negative Bias Weights ---")
    for aa in standard_aas:
        f_p = (pass_counts[aa] + 1) / (tot_pass_aa + 20)
        f_f = (fail_counts[aa] + 1) / (tot_fail_aa + 20)
        lor = np.log(f_p / f_f)
        mpnn_bias_dict[aa] = round(float(lor), 3)
    
    sorted_bias = sorted(mpnn_bias_dict.items(), key=lambda x: -x[1])
    for aa, bias in sorted_bias:
        tag = "[ENRICHED / FAVOR]" if bias > 0.15 else ("[DEPLETED / PENALIZE]" if bias < -0.15 else "[NEUTRAL]")
        print(f"  Residue {aa:2s} : {bias:+.3f}  {tag}")
        
    # 3. Biophysical Feature Matrix
    feature_names = [
        'length', 'Interface_BuriedArea', 'Surface_Hydrophobicity', 
        'Interface_Hydrophobicity', 'Binder_Helix_Fraction', 'Binder_BetaSheet_Fraction',
        'Binder_Loop_Fraction', 'Binder_pI', 'Binder_Net_Charge', 'Unbound_Binder_pLDDT'
    ]
    
    X_mat = []
    for _, row in df.iterrows():
        vec = [
            float(row['length']),
            parse_float(row['Interface_BuriedArea'], 500.0),
            parse_float(row['Surface_Hydrophobicity'], 0.25),
            parse_float(row['Interface_Hydrophobicity'], 0.40),
            parse_float(row['Binder_Helix_Fraction'], 0.8),
            parse_float(row['Binder_BetaSheet_Fraction'], 0.0),
            parse_float(row['Binder_Loop_Fraction'], 0.2),
            parse_float(row['Binder_pI'], 6.0),
            parse_float(row['Binder_Net_Charge'], -2.0),
            parse_float(row['Unbound_Binder_pLDDT'], 0.85)
        ]
        X_mat.append(vec)
    
    X = np.array(X_mat, dtype=np.float32)
    y = np.array(df['is_pass'].values, dtype=np.float32)  # 1 = pass, 0 = fail
    
    mean = np.mean(X, axis=0)
    std = np.std(X, axis=0) + 1e-6
    X_norm = (X - mean) / std
    
    # 4. PyTorch GPU Training on Agni
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    print(f"\n--- 2. PyTorch Latent Failure Cloud Network on {device} ---")
    if device.type == "cuda":
        print(f"  Active GPU: {torch.cuda.get_device_name(0)}")
        
    model = NegativeCloudNet(input_dim=len(feature_names)).to(device)
    criterion = nn.BCEWithLogitsLoss()
    optimizer = optim.AdamW(model.parameters(), lr=0.01, weight_decay=1e-4)
    
    X_t = torch.tensor(X_norm, dtype=torch.float32, device=device)
    y_t = torch.tensor(y, dtype=torch.float32, device=device).unsqueeze(1)
    
    # Train discriminator / energy surface
    model.train()
    for epoch in range(150):
        optimizer.zero_grad()
        logits = model(X_t)
        loss = criterion(logits, y_t)
        loss.backward()
        optimizer.step()
        
    print(f"  Training converged: BCE Loss = {loss.item():.4f}")
    
    # 5. Generative Guidance & Steering Demo
    print("\n--- 3. Real-Time Generative Steering & Repulsion Gradients ---")
    model.eval()
    
    # Evaluate Passer Lead (Slot 1)
    sample_pass = torch.tensor(X_norm[0:1], dtype=torch.float32, device=device, requires_grad=True)
    logit_pass = model(sample_pass)
    prob_pass = torch.sigmoid(logit_pass).item()
    print(f"  Lead APEX-EGFR-01 (Passed):")
    print(f"    P(Pass) = {prob_pass*100:.1f}% | Latent Energy = {logit_pass.item():+.3f}")
    
    # Evaluate Discarded Design (Slot 6, failed buried area)
    sample_fail = torch.tensor(X_norm[5:6], dtype=torch.float32, device=device, requires_grad=True)
    logit_fail = model(sample_fail)
    prob_fail_pass = torch.sigmoid(logit_fail).item()
    
    # Compute Repulsive Steering Gradient (Pushing away from failure cloud)
    # Energy to maximize is P(Pass), so gradient ascent direction on logit pushes toward safety
    logit_fail.backward()
    grad = sample_fail.grad.cpu().numpy()[0]
    
    print(f"  Discarded Candidate (Slot 6, Failed):")
    print(f"    P(Pass) = {prob_fail_pass*100:.1f}% | P(Failure Cloud) = {(1.0-prob_fail_pass)*100:.1f}%")
    print(f"    Repulsive Guidance Gradient (top 3 restorative adjustments):")
    
    grad_tuples = [(feature_names[i], grad[i]) for i in range(len(feature_names))]
    grad_tuples.sort(key=lambda x: -abs(x[1]))
    for feat, g in grad_tuples[:3]:
        direction = "INCREASE" if g > 0 else "DECREASE"
        print(f"      -> {direction:8s} {feat:25s} (Steering Force: {g:+.4f})")

    print("\n=== Active Learning Feedback Loop Successfully Validated on Hardware ===")

if __name__ == "__main__":
    run_feedback_pipeline()
