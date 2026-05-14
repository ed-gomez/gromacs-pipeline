# GROMACS Protein MD Simulation Script
A bash script that will prepare a pdb file for a MD run and handle minimization, equilibration, and production runs automatically.

---

## Quick Start

```bash
# Download 1AKI and run with all defaults
./run_sim.sh 1AKI

# Run a local PDB file
./run_sim.sh myprotein.pdb

# Longer production run at higher temperature
./run_sim.sh --md-ns 100 --temp 310 2LZM.pdb

# Resume an interrupted run from NVT
./run_sim.sh --resume nvt 1AKI
```

---

## Pipeline Overview

```
Raw PDB
  │
  ▼ grep ATOM (remove waters/ligands)
Clean PDB
  │
  ▼ gmx pdb2gmx   (-ff oplsaa -water spce)
  │   Assigns force field parameters, protonation states, hydrogen atoms
  │   → topol.top, posre.itp, *_processed.gro
  │
  ▼ gmx editconf   (-d 1.0 -bt cubic)
  │   Defines a periodic simulation box with ≥1 nm padding around protein
  │   → *_newbox.gro
  │
  ▼ gmx solvate    (-cs spc216.gro)
  │   Fills box with explicit water molecules
  │   → *_solv.gro  (topol.top updated)
  │
  ▼ gmx grompp + gmx genion   (-neutral)
  │   Adds Na⁺/Cl⁻ ions to neutralize net charge of protein
  │   → *_solv_ions.gro  (topol.top updated)
  │
  ▼ Energy Minimization   (steepest descent, ≤50,000 steps)
  │   Relaxes bad contacts introduced by solvation/ion placement
  │   → em.gro, em.edr, em.log
  │
  ▼ NVT Equilibration   (100 ps, 300 K, position-restrained)
  │   Heats the system to target temperature with protein held in place
  │   → nvt.gro, nvt.cpt, nvt.edr
  │
  ▼ NPT Equilibration   (100 ps, 300 K / 1 bar, position-restrained)
  │   Equilibrates pressure and density with protein held in place
  │   → npt.gro, npt.cpt, npt.edr
  │
  ▼ Production MD   (1 ns, no restraints)
      Free dynamics — the simulation data you actually analyze
      → md.xtc, md.gro, md.edr
```

---

## Key Parameters

### Force field (`-ff`)
The force field defines atomic charges, bond lengths, angles, and torsions.

| Flag value | Force field | Best for |
|---|---|---|
| `oplsaa` | OPLS-AA (default) | small organic molecules, good general choice |
| `amber99sb-ildn` | AMBER99SB-ILDN | proteins, widely benchmarked |
| `charmm36` | CHARMM36 | proteins, membranes |

Run `gmx pdb2gmx` with no flags to see all available force fields.

### Water model (`-w`)
Must be compatible with the chosen force field.

| Flag | Model | Notes |
|---|---|---|
| `spce` | SPC/E (default) | Good general choice, compatible with OPLS-AA |
| `tip3p` | TIP3P | Required with CHARMM36/AMBER |
| `tip4p` | TIP4P | More accurate water properties |

### Box type (`-b`)
| Flag | Shape | Volume efficiency | Notes |
|---|---|---|---|
| `cubic` | Cube | ~100% | Simplest, largest box |
| `dodecahedron` | Rhombic dodecahedron | ~71% | Good default for globular proteins |
| `octahedron` | Truncated octahedron | ~77% | Fewer water molecules than cubic |

Using `dodecahedron` instead of `cubic` with the same padding distance reduces water molecules by ~30%, making simulations ~30% faster.

### Simulation lengths
The defaults (NVT=0.1ns, NPT=0.1ns, MD=1ns) are fine for tutorial/testing.

| Use case | NVT | NPT | MD |
|---|---|---|---|
| Tutorial / testing | 0.1 ns | 0.1 ns | 1 ns |
| Structural analysis | 0.1 ns | 0.1 ns | 10–100 ns |
| Free energy / binding | 0.5 ns | 0.5 ns | 100–1000 ns |

### Temperature (`-T`)
Physiological temperature is 310 K (37°C). The tutorial uses 300 K (27°C) — common for benchmarking.

---

## How the MDP Files Work

Each simulation phase has a parameter file (`.mdp`). The script generates these automatically based on your options. Key settings:

### minim.mdp — energy minimization
- `integrator = steep` — steepest-descent algorithm (not MD)
- `emtol = 1000.0` — stop when max force < 1000 kJ/mol/nm
- No temperature or pressure coupling (not a dynamics run)

### nvt.mdp — NVT equilibration
- `define = -DPOSRES` — activates position restraints on protein heavy atoms (from `posre.itp`)
- `tcoupl = V-rescale` — velocity rescaling thermostat
- `pcoupl = no` — no pressure coupling in NVT
- `gen_vel = yes` — assigns velocities from Maxwell-Boltzmann distribution

### npt.mdp — NPT equilibration
- `define = -DPOSRES` — still position-restrained
- `pcoupl = Parrinello-Rahman` — barostat for pressure control
- `continuation = yes` — reads velocities from NVT checkpoint
- `gen_vel = no` — velocities already set

### md.mdp — production MD
- No position restraints (protein moves freely)
- `nstxout = 0` — suppresses bulky `.trr` file; only compressed `.xtc` is written
- `nstxout-compressed = 5000` — saves frame every 10 ps

---

## Analyzing Your Results

All commands run from inside the protein's output directory (e.g., `1aki/`):

### Verify energy minimization
```bash
echo "Potential" | gmx energy -f em.edr -o potential.xvg
xmgrace potential.xvg   # or: python plot.py potential.xvg
```
The potential should decrease monotonically and plateau.

### Check NVT temperature stabilization
```bash
echo "Temperature" | gmx energy -f nvt.edr -o temperature.xvg
```
Temperature should converge to your target (e.g., 300 K) within ~20 ps.

### Check NPT pressure and density
```bash
echo "Pressure" | gmx energy -f npt.edr -o pressure.xvg
echo "Density"  | gmx energy -f npt.edr -o density.xvg
```
Density for TIP3P/SPC-E water should converge to ~1000 kg/m³.

### RMSD — structural drift during production
```bash
echo "Backbone
Backbone" | gmx rms -s md.tpr -f md.xtc -o rmsd.xvg -tu ns
```
Low RMSD (<0.2 nm for a rigid protein) means the structure is stable.

### RMSF — per-residue flexibility
```bash
echo "C-alpha" | gmx rmsf -s md.tpr -f md.xtc -o rmsf.xvg -res
```
High RMSF residues are flexible loops or termini.

### Radius of gyration — compactness
```bash
echo "Protein" | gmx gyrate -s md.tpr -f md.xtc -o gyrate.xvg
```

### Prepare trajectory for visualization (VMD/PyMOL)
```bash
# Remove periodic boundary artifacts, center protein
echo "Protein
System" | gmx trjconv -s md.tpr -f md.xtc -o md_center.xtc -pbc mol -center

# Or extract just the first frame as a PDB
echo "System" | gmx trjconv -s md.tpr -f md.xtc -o frame0.pdb -dump 0
```

---

## Common Issues

### `Fatal error: atom XX not found in residue YYY`
The PDB has a non-standard residue or unusual atom naming. Options:
- Use `grep "^ATOM"` (already done by script) to remove ligands, then handle them separately
- Try a different force field
- Manually rename atoms to match the force field's expected names

### `WARNING: 1 bonds with force constants > 20000`
Usually harmless for the tutorial. Add `-maxwarn 1` to the `grompp` call if needed.

### `LINCS warning` during MD
Hydrogen bond constraints failed. This often means the system wasn't equilibrated enough. Try:
- Longer NVT/NPT equilibration (`--nvt-ns 0.5 --npt-ns 0.5`)
- Check for unrealistic clashes after genion step (view `*_solv_ions.gro` in VMD)

### `genion: group SOL not found`
The script pipes `"SOL"` to genion. If your topology uses a different name, edit the genion call in the script to pipe the correct group number (run `gmx genion -s ions.tpr` interactively first to see the group list).

### Simulation is very slow
Pass threading/GPU arguments:
```bash
./run_sim.sh --mdrun-args "-ntmpi 1 -ntomp 8 -gpu_id 0" 1AKI
```

---

## File Reference

| File | Created by | Contents |
|---|---|---|
| `*_clean.pdb` | script | PDB with only ATOM records |
| `*_processed.gro` | pdb2gmx | Protein structure with H atoms |
| `topol.top` | pdb2gmx | Full system topology |
| `posre.itp` | pdb2gmx | Position restraint parameters |
| `*_newbox.gro` | editconf | Protein in periodic box |
| `*_solv.gro` | solvate | Protein + water |
| `*_solv_ions.gro` | genion | Protein + water + ions |
| `em.gro/edr/log` | mdrun (EM) | Minimized structure + energy |
| `nvt.gro/cpt/edr` | mdrun (NVT) | NVT-equilibrated structure |
| `npt.gro/cpt/edr` | mdrun (NPT) | NPT-equilibrated structure |
| `md.xtc` | mdrun (MD) | Production trajectory (compressed) |
| `md.gro` | mdrun (MD) | Final structure |
| `md.edr` | mdrun (MD) | Energy/thermodynamic data |

---

## Reference

- Tutorial: http://www.mdtutorials.com/gmx/lysozyme/
- GROMACS manual: https://manual.gromacs.org/
- GROMACS 2024 release notes: https://manual.gromacs.org/2024/release-notes/
