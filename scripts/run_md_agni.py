import os
import sys
import openmm
from openmm import app, unit
from pdbfixer import PDBFixer

def run_simulation(pdb_file, ph_val, out_prefix, steps=5000000): # 10 ns = 5,000,000 steps of 2 fs
    print(f"=== Starting MD Simulation: {out_prefix} at pH {ph_val} ===", flush=True)
    os.makedirs("/home/agni/ApexEGFR/md_simulations/outputs", exist_ok=True)
    
    # 1. PDB Fixer
    print(f"1. Fixing PDB topology and protonation at pH {ph_val}...", flush=True)
    fixer = PDBFixer(filename=pdb_file)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH=ph_val)
    print(f"   Atoms after PDBFixer: {fixer.topology.getNumAtoms()}", flush=True)
    
    # 2. Forcefield and Solvation
    print("2. Building TIP3P water box (1.0 nm padding, 150 mM NaCl)...", flush=True)
    ff = app.ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
    modeller = app.Modeller(fixer.topology, fixer.positions)
    modeller.addSolvent(ff, model='tip3p', padding=1.0*unit.nanometers, ionicStrength=0.15*unit.molar)
    print(f"   Total solvated atoms: {modeller.topology.getNumAtoms()}", flush=True)
    
    # 3. System creation
    print("3. Creating OpenMM System (PME, 1.0 nm cutoff, HBonds constraints)...", flush=True)
    system = ff.createSystem(
        modeller.topology,
        nonbondedMethod=app.PME,
        nonbondedCutoff=1.0*unit.nanometers,
        constraints=app.HBonds
    )
    system.addForce(openmm.MonteCarloBarostat(1.0*unit.atmospheres, 310.15*unit.kelvin, 25))
    
    # 4. Integrator and Platform
    integrator = openmm.LangevinMiddleIntegrator(
        310.15*unit.kelvin,
        1.0/unit.picoseconds,
        2.0*unit.femtoseconds
    )
    platform = openmm.Platform.getPlatformByName("OpenCL")
    properties = {'OpenCLPrecision': 'mixed'}
    
    simulation = app.Simulation(modeller.topology, system, integrator, platform, properties)
    simulation.context.setPositions(modeller.positions)
    
    # 5. Energy Minimization
    print("4. Energy Minimization (Steepest Descent)...", flush=True)
    simulation.minimizeEnergy(maxIterations=1000)
    print("   Minimization complete!", flush=True)
    
    # 6. Equilibrating NVT / NPT (50 ps = 25,000 steps)
    print("5. Equilibration (50 ps at 310.15 K)...", flush=True)
    simulation.step(25000)
    print("   Equilibration complete!", flush=True)
    
    # 7. Production Run with Reporters
    dcd_path = f"/home/agni/ApexEGFR/md_simulations/outputs/{out_prefix}.dcd"
    log_path = f"/home/agni/ApexEGFR/md_simulations/outputs/{out_prefix}.csv"
    pdb_equil = f"/home/agni/ApexEGFR/md_simulations/outputs/{out_prefix}_equilibrated.pdb"
    
    with open(pdb_equil, 'w') as f:
        app.PDBFile.writeFile(simulation.topology, simulation.context.getState(getPositions=True).getPositions(), f)
        
    print(f"6. Production Run ({steps} steps = {steps*2/1e6:.1f} ns) -> logging to {log_path}", flush=True)
    simulation.reporters.append(app.DCDReporter(dcd_path, 10000)) # Save frame every 20 ps
    simulation.reporters.append(app.StateDataReporter(
        log_path, 5000, step=True, potentialEnergy=True, temperature=True, speed=True
    ))
    simulation.reporters.append(app.StateDataReporter(
        sys.stdout, 5000, step=True, potentialEnergy=True, temperature=True, speed=True
    ))
    
    simulation.step(steps)
    print(f"=== Simulation {out_prefix} Finished Successfully! ===", flush=True)

if __name__ == '__main__':
    pdb = '/home/agni/ApexEGFR/md_simulations/apex_egfr_02_hEGFR.pdb'
    mode = sys.argv[1] if len(sys.argv) > 1 else 'ph65'
    
    if mode == 'ph65':
        run_simulation(pdb, 6.5, 'apex_02_tumor_pH65', steps=5000000)
    elif mode == 'ph74':
        run_simulation(pdb, 7.4, 'apex_02_physio_pH74', steps=5000000)
