#!/usr/bin/env bash
# run_sim.sh — Automate a standard GROMACS protein MD simulation
#
# Pipeline: clean PDB → pdb2gmx → editconf → solvate → genion
#           → energy minimization → NVT equilibration → NPT equilibration → production MD
#
# Usage: ./run_sim.sh [OPTIONS] <protein.pdb | PDBID>
# See --help for full options.

set -euo pipefail

# ── defaults ──────────────────────────────────────────────────────────────────
FF="oplsaa"
WATER="spce"
BOX_DIST="1.0"       # nm from protein to box edge
BOX_TYPE="cubic"
TEMP="300"           # K
PRESSURE="1.0"       # bar
NVT_NS="0.1"         # 100 ps
NPT_NS="0.1"         # 100 ps
MD_NS="1.0"          # 1 ns
DT="0.002"           # ps timestep (2 fs)
PNAME="NA"
NNAME="CL"
GMX="gmx"
RESUME_FROM=""
EXTRA_MDRUN=""

# ── helpers ───────────────────────────────────────────────────────────────────
die()  { echo "ERROR: $*" >&2; exit 1; }
info() { echo "[$(date '+%H:%M:%S')] >>> $*"; }
ok()   { echo "[$(date '+%H:%M:%S')]  ✓  $*"; }

ns_to_steps() {
    # Convert nanoseconds to MD steps using global $DT (ps)
    printf "%.0f" "$(echo "$1 * 1000 / $DT" | bc -l)"
}

# Map water model to solvate template file
water_template() {
    case "$1" in
        tip4p|tip4pew) echo "tip4p.gro" ;;
        *)             echo "spc216.gro" ;;
    esac
}

usage() {
    cat <<'EOF'
Usage: run_sim.sh [OPTIONS] <protein.pdb | PDBID>

Automates the standard GROMACS protein MD simulation pipeline:
  PDB prep → topology → solvation → ions → EM → NVT → NPT → MD

Arguments:
  protein.pdb   Path to a PDB file
  PDBID         4-character RCSB PDB ID (e.g. 1AKI) — downloaded automatically

Options:
  -ff, --forcefield FF   Force field name passed to pdb2gmx (default: oplsaa)
  -w,  --water MODEL     Water model: spce, tip3p, tip4p (default: spce)
  -d,  --box-dist NM     Min distance protein-to-box edge in nm (default: 1.0)
  -b,  --box-type TYPE   Box type: cubic, dodecahedron, octahedron (default: cubic)
  -T,  --temp K          Simulation temperature in K (default: 300)
  -P,  --pressure BAR    Reference pressure in bar (default: 1.0)
       --nvt-ns NS       NVT equilibration length in ns (default: 0.1)
       --npt-ns NS       NPT equilibration length in ns (default: 0.1)
       --md-ns  NS       Production MD length in ns (default: 1.0)
       --gmx PATH        Path to gmx executable (default: gmx)
       --resume STEP     Skip steps before STEP and resume from there.
                         STEP: pdb2gmx | box | solvate | ions | em | nvt | npt | md
       --mdrun-args "…"  Extra arguments forwarded to every mdrun call
                         Example: --mdrun-args "-ntmpi 1 -ntomp 8 -gpu_id 0"
  -h,  --help            Show this help

Examples:
  # Run 1AKI from scratch, downloading the PDB automatically
  ./run_sim.sh 1AKI

  # Run a local PDB with a longer production and higher temperature
  ./run_sim.sh --md-ns 10 --temp 310 myprotein.pdb

  # Resume an interrupted run from the NPT step
  ./run_sim.sh --resume npt 1AKI

  # Use GPU and multiple threads
  ./run_sim.sh --mdrun-args "-ntmpi 1 -ntomp 8 -gpu_id 0" 1AKI

Output:
  A subdirectory named after the protein (e.g. ./1aki/) containing all outputs.
EOF
    exit 0
}

# ── argument parsing ──────────────────────────────────────────────────────────
[[ $# -eq 0 ]] && usage

INPUT=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        -ff|--forcefield)  FF="$2";            shift 2 ;;
        -w|--water)        WATER="$2";         shift 2 ;;
        -d|--box-dist)     BOX_DIST="$2";      shift 2 ;;
        -b|--box-type)     BOX_TYPE="$2";      shift 2 ;;
        -T|--temp)         TEMP="$2";          shift 2 ;;
        -P|--pressure)     PRESSURE="$2";      shift 2 ;;
        --nvt-ns)          NVT_NS="$2";        shift 2 ;;
        --npt-ns)          NPT_NS="$2";        shift 2 ;;
        --md-ns)           MD_NS="$2";         shift 2 ;;
        --gmx)             GMX="$2";           shift 2 ;;
        --resume)          RESUME_FROM="$2";   shift 2 ;;
        --mdrun-args)      EXTRA_MDRUN="$2";   shift 2 ;;
        -h|--help)         usage ;;
        -*)                die "Unknown option: $1. Run with --help." ;;
        *)                 INPUT="$1";         shift ;;
    esac
done

[[ -z "$INPUT" ]] && die "No PDB file or PDB ID provided. Run with --help."

# Validate gmx
command -v "$GMX" &>/dev/null || die "'$GMX' not found. Install GROMACS or set --gmx /path/to/gmx"

# ── derive protein name and working directory ─────────────────────────────────
if [[ "$INPUT" =~ ^[A-Za-z0-9]{4}$ && ! -f "$INPUT" ]]; then
    PDBID=$(echo "$INPUT" | tr '[:lower:]' '[:upper:]')
    NAME=$(echo "$INPUT"  | tr '[:upper:]' '[:lower:]')
    DOWNLOAD_PDB=true
    ABS_INPUT=""
elif [[ -f "$INPUT" ]]; then
    PDBID=""
    NAME=$(basename "$INPUT" .pdb)
    NAME=$(basename "$NAME" .PDB)
    NAME=$(echo "$NAME" | tr '[:upper:]' '[:lower:]')
    DOWNLOAD_PDB=false
    ABS_INPUT="$(cd "$(dirname "$INPUT")"; pwd)/$(basename "$INPUT")"
else
    die "'$INPUT' is not an existing file and doesn't look like a 4-char PDB ID."
fi

WORKDIR="$(pwd)/${NAME}"
mkdir -p "$WORKDIR"
cd "$WORKDIR"

LOGFILE="${WORKDIR}/${NAME}_sim.log"
# Tee all output to log file
exec > >(tee -a "$LOGFILE") 2>&1

info "Starting GROMACS simulation for: ${NAME}"
info "Working directory: ${WORKDIR}"
info "Force field: ${FF} | Water: ${WATER} | Box: ${BOX_TYPE} d=${BOX_DIST}nm"
info "Temperature: ${TEMP}K | Pressure: ${PRESSURE}bar"
info "NVT: ${NVT_NS}ns | NPT: ${NPT_NS}ns | MD: ${MD_NS}ns"
[[ -n "$RESUME_FROM" ]] && info "Resuming from step: ${RESUME_FROM}"

# ── step ordering for resume logic ────────────────────────────────────────────
STEPS=(pdb2gmx box solvate ions em nvt npt md)

step_index() {
    local target="$1"
    local i=0
    for s in "${STEPS[@]}"; do
        [[ "$s" == "$target" ]] && echo "$i" && return
        ((i++))
    done
    die "Unknown step '$target'. Valid steps: ${STEPS[*]}"
}

RESUME_IDX=0
[[ -n "$RESUME_FROM" ]] && RESUME_IDX=$(step_index "$RESUME_FROM")

should_run() {
    local idx
    idx=$(step_index "$1")
    [[ "$idx" -ge "$RESUME_IDX" ]]
}

# ── MDP file writers ──────────────────────────────────────────────────────────
write_ions_mdp() {
    cat > ions.mdp <<EOF
; ions.mdp — single-point calculation to generate ions.tpr
integrator      = steep
emtol           = 1000.0
emstep          = 0.01
nsteps          = 50000
nstlist         = 1
cutoff-scheme   = Verlet
ns_type         = grid
coulombtype     = cutoff
rcoulomb        = 1.0
rvdw            = 1.0
pbc             = xyz
EOF
}

write_minim_mdp() {
    cat > minim.mdp <<EOF
; minim.mdp — steepest-descent energy minimization
integrator      = steep
emtol           = 1000.0
emstep          = 0.01
nsteps          = 50000
nstlist         = 1
cutoff-scheme   = Verlet
ns_type         = grid
coulombtype     = PME
rcoulomb        = 1.0
rvdw            = 1.0
pbc             = xyz
EOF
}

write_nvt_mdp() {
    local steps
    steps=$(ns_to_steps "$NVT_NS")
    cat > nvt.mdp <<EOF
; nvt.mdp — NVT equilibration (${NVT_NS} ns, position-restrained)
define                  = -DPOSRES
integrator              = md
nsteps                  = ${steps}
dt                      = ${DT}
nstxout                 = 500
nstvout                 = 500
nstenergy               = 500
nstlog                  = 500
continuation            = no
constraint_algorithm    = lincs
constraints             = h-bonds
lincs_iter              = 1
lincs_order             = 4
cutoff-scheme           = Verlet
ns_type                 = grid
nstlist                 = 10
rcoulomb                = 1.0
rvdw                    = 1.0
DispCorr                = EnerPres
coulombtype             = PME
pme_order               = 4
fourierspacing          = 0.16
tcoupl                  = V-rescale
tc-grps                 = Protein Non-Protein
tau_t                   = 0.1     0.1
ref_t                   = ${TEMP}     ${TEMP}
pcoupl                  = no
pbc                     = xyz
gen_vel                 = yes
gen_temp                = ${TEMP}
gen_seed                = -1
EOF
}

write_npt_mdp() {
    local steps
    steps=$(ns_to_steps "$NPT_NS")
    cat > npt.mdp <<EOF
; npt.mdp — NPT equilibration (${NPT_NS} ns, position-restrained)
define                  = -DPOSRES
integrator              = md
nsteps                  = ${steps}
dt                      = ${DT}
nstxout                 = 500
nstvout                 = 500
nstenergy               = 500
nstlog                  = 500
continuation            = yes
constraint_algorithm    = lincs
constraints             = h-bonds
lincs_iter              = 1
lincs_order             = 4
cutoff-scheme           = Verlet
ns_type                 = grid
nstlist                 = 10
rcoulomb                = 1.0
rvdw                    = 1.0
DispCorr                = EnerPres
coulombtype             = PME
pme_order               = 4
fourierspacing          = 0.16
tcoupl                  = V-rescale
tc-grps                 = Protein Non-Protein
tau_t                   = 0.1     0.1
ref_t                   = ${TEMP}     ${TEMP}
pcoupl                  = Parrinello-Rahman
pcoupltype              = isotropic
tau_p                   = 2.0
ref_p                   = ${PRESSURE}
compressibility         = 4.5e-5
refcoord_scaling        = com
pbc                     = xyz
gen_vel                 = no
EOF
}

write_md_mdp() {
    local steps
    steps=$(ns_to_steps "$MD_NS")
    cat > md.mdp <<EOF
; md.mdp — production MD (${MD_NS} ns)
integrator              = md
nsteps                  = ${steps}
dt                      = ${DT}
nstxout                 = 0
nstvout                 = 0
nstfout                 = 0
nstenergy               = 5000
nstlog                  = 5000
nstxout-compressed      = 5000
compressed-x-grps       = System
continuation            = yes
constraint_algorithm    = lincs
constraints             = h-bonds
lincs_iter              = 1
lincs_order             = 4
cutoff-scheme           = Verlet
ns_type                 = grid
nstlist                 = 10
rcoulomb                = 1.0
rvdw                    = 1.0
DispCorr                = EnerPres
coulombtype             = PME
pme_order               = 4
fourierspacing          = 0.16
tcoupl                  = V-rescale
tc-grps                 = Protein Non-Protein
tau_t                   = 0.1     0.1
ref_t                   = ${TEMP}     ${TEMP}
pcoupl                  = Parrinello-Rahman
pcoupltype              = isotropic
tau_p                   = 2.0
ref_p                   = ${PRESSURE}
compressibility         = 4.5e-5
pbc                     = xyz
gen_vel                 = no
EOF
}

# ── STEP 0: obtain PDB ────────────────────────────────────────────────────────
if $DOWNLOAD_PDB; then
    RAW_PDB="${NAME}_raw.pdb"
    if [[ ! -f "$RAW_PDB" ]]; then
        info "Downloading ${PDBID} from RCSB..."
        curl -fsSL "https://files.rcsb.org/download/${PDBID}.pdb" -o "$RAW_PDB" \
            || die "Download failed. Check the PDB ID and your internet connection."
        ok "Downloaded ${PDBID}.pdb"
    else
        info "Found existing ${RAW_PDB}, skipping download."
    fi
    INPUT_PDB="$RAW_PDB"
else
    RAW_PDB="${NAME}_raw.pdb"
    if [[ ! -f "$RAW_PDB" ]]; then
        cp "$ABS_INPUT" "$RAW_PDB" || die "Cannot copy PDB file: $ABS_INPUT"
    fi
    INPUT_PDB="$RAW_PDB"
fi

# Clean PDB: keep ATOM records only (removes HOH/crystallographic waters and HETATM ligands)
CLEAN_PDB="${NAME}_clean.pdb"
if [[ ! -f "$CLEAN_PDB" ]]; then
    info "Cleaning PDB (removing HETATM records)..."
    grep "^ATOM" "$INPUT_PDB" > "$CLEAN_PDB"
    ATOM_COUNT=$(wc -l < "$CLEAN_PDB")
    ok "Cleaned PDB: ${ATOM_COUNT} ATOM lines → ${CLEAN_PDB}"
else
    info "Found ${CLEAN_PDB}, skipping clean step."
fi

[[ ! -s "$CLEAN_PDB" ]] && die "${CLEAN_PDB} is empty — check that the PDB file has ATOM records."

# ── STEP 1: pdb2gmx ───────────────────────────────────────────────────────────
if should_run pdb2gmx; then
    info "Step 1/8: pdb2gmx — generating topology and processed structure"
    $GMX pdb2gmx \
        -f "$CLEAN_PDB" \
        -o "${NAME}_processed.gro" \
        -p topol.top \
        -i posre.itp \
        -ff "$FF" \
        -water "$WATER" \
        -ignh
    ok "Topology and processed structure written."
else
    info "Skipping pdb2gmx (resume mode)."
fi

# ── STEP 2: editconf — define simulation box ──────────────────────────────────
if should_run box; then
    info "Step 2/8: editconf — defining simulation box (${BOX_TYPE}, d=${BOX_DIST}nm)"
    $GMX editconf \
        -f "${NAME}_processed.gro" \
        -o "${NAME}_newbox.gro" \
        -c \
        -d "$BOX_DIST" \
        -bt "$BOX_TYPE"
    ok "Simulation box defined."
else
    info "Skipping editconf (resume mode)."
fi

# ── STEP 3: solvate ───────────────────────────────────────────────────────────
if should_run solvate; then
    SOLVENT_TEMPLATE=$(water_template "$WATER")
    info "Step 3/8: solvate — filling box with water (${WATER} / ${SOLVENT_TEMPLATE})"
    $GMX solvate \
        -cp "${NAME}_newbox.gro" \
        -cs "$SOLVENT_TEMPLATE" \
        -o "${NAME}_solv.gro" \
        -p topol.top
    ok "System solvated."
else
    info "Skipping solvate (resume mode)."
fi

# ── STEP 4: genion — add ions to neutralize ───────────────────────────────────
if should_run ions; then
    info "Step 4/8: genion — adding counter-ions to neutralize system"
    write_ions_mdp
    $GMX grompp \
        -f ions.mdp \
        -c "${NAME}_solv.gro" \
        -p topol.top \
        -o ions.tpr \
        -maxwarn 2
    echo "SOL" | $GMX genion \
        -s ions.tpr \
        -o "${NAME}_solv_ions.gro" \
        -p topol.top \
        -pname "$PNAME" \
        -nname "$NNAME" \
        -neutral
    ok "Ions added. System is now charge-neutral."
else
    info "Skipping genion (resume mode)."
fi

# ── STEP 5: energy minimization ───────────────────────────────────────────────
if should_run em; then
    info "Step 5/8: energy minimization (steepest descent, max 50000 steps)"
    write_minim_mdp
    $GMX grompp \
        -f minim.mdp \
        -c "${NAME}_solv_ions.gro" \
        -p topol.top \
        -o em.tpr
    # shellcheck disable=SC2086
    $GMX mdrun -v -deffnm em $EXTRA_MDRUN

    # Check convergence
    CONVERGED=$(grep "converged to Fmax" em.log 2>/dev/null | tail -1 || true)
    if [[ -z "$CONVERGED" ]]; then
        echo "WARNING: Energy minimization may not have converged. Check em.log." >&2
    else
        ok "Energy minimization converged: ${CONVERGED}"
    fi
else
    info "Skipping energy minimization (resume mode)."
fi

# ── STEP 6: NVT equilibration ─────────────────────────────────────────────────
if should_run nvt; then
    info "Step 6/8: NVT equilibration (${NVT_NS} ns at ${TEMP}K, position-restrained)"
    write_nvt_mdp
    $GMX grompp \
        -f nvt.mdp \
        -c em.gro \
        -r em.gro \
        -p topol.top \
        -o nvt.tpr
    # shellcheck disable=SC2086
    $GMX mdrun -deffnm nvt $EXTRA_MDRUN
    ok "NVT equilibration complete."
else
    info "Skipping NVT (resume mode)."
fi

# ── STEP 7: NPT equilibration ─────────────────────────────────────────────────
if should_run npt; then
    info "Step 7/8: NPT equilibration (${NPT_NS} ns at ${TEMP}K / ${PRESSURE}bar, position-restrained)"
    write_npt_mdp
    $GMX grompp \
        -f npt.mdp \
        -c nvt.gro \
        -r nvt.gro \
        -t nvt.cpt \
        -p topol.top \
        -o npt.tpr
    # shellcheck disable=SC2086
    $GMX mdrun -deffnm npt $EXTRA_MDRUN
    ok "NPT equilibration complete."
else
    info "Skipping NPT (resume mode)."
fi

# ── STEP 8: production MD ─────────────────────────────────────────────────────
if should_run md; then
    info "Step 8/8: production MD (${MD_NS} ns)"
    write_md_mdp
    $GMX grompp \
        -f md.mdp \
        -c npt.gro \
        -t npt.cpt \
        -p topol.top \
        -o md.tpr
    # shellcheck disable=SC2086
    $GMX mdrun -deffnm md $EXTRA_MDRUN
    ok "Production MD complete."
else
    info "Skipping production MD (resume mode)."
fi

# ── done ──────────────────────────────────────────────────────────────────────
cat <<EOF

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  Simulation complete: ${NAME}
  Output directory:    ${WORKDIR}/
  Trajectory:          md.xtc
  Final structure:     md.gro
  Full log:            ${LOGFILE}
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Quick analysis commands (run from ${WORKDIR}/):

  # Potential energy during EM
  echo "Potential" | gmx energy -f em.edr -o potential.xvg

  # Temperature during NVT
  echo "Temperature" | gmx energy -f nvt.edr -o temperature.xvg

  # Pressure during NPT
  echo "Pressure" | gmx energy -f npt.edr -o pressure.xvg

  # Density during NPT
  echo "Density" | gmx energy -f npt.edr -o density.xvg

  # RMSD vs initial structure (production MD)
  echo "Backbone Backbone" | gmx rms -s md.tpr -f md.xtc -o rmsd.xvg -tu ns

  # Radius of gyration
  echo "Protein" | gmx gyrate -s md.tpr -f md.xtc -o gyrate.xvg

  # Strip periodic boundary artifacts for visualization
  echo "Protein" | gmx trjconv -s md.tpr -f md.xtc -o md_noPBC.xtc -pbc mol -center

EOF
