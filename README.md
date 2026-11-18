# Center_Trajectory

A Python 3 command-line script that reconstructs molecules across periodic boundaries, centers a user-selected atom group in the simulation box, wraps molecules back into the box, and writes a trajectory using a possibly different atom selection.

## Requirements

- Python 3
- MDAnalysis 2.4 or newer
- A full-system topology readable by MDAnalysis that includes bonds, such as a GROMACS `.tpr`
- A full-system GROMACS `.xtc` trajectory containing the same atoms as the structure

A `.gro` file alone does not contain the bond information needed to reconstruct whole molecules. Use a topology-bearing input such as `.tpr`.

Create the Conda environment from the included YAML file:

```bash
conda env create -f environment.yaml
conda activate center-protein
```

Or install MDAnalysis into an existing Python 3 environment:

```bash
python3 -m pip install -r requirements.txt
```

## What it does

For each trajectory frame, the script:

1. Uses MDAnalysis `unwrap` to reconstruct molecules from the topology's bonds.
2. Centers the group specified by `--selection` in the simulation box.
3. Uses MDAnalysis `wrap` to place whole molecules in the periodic box.
4. Writes atoms matched by `--output-selection` to the output XTC.

The default centering selection is `protein`. If `--output-selection` is omitted, the script writes the same group used for centering.

## Usage

```bash
python3 center.py \
  --structure system.tpr \
  --trajectory trajectory.xtc \
  --selection protein \
  --output-selection protein \
  --output protein_centered.xtc \
  --output-gro protein_centered_final.gro
```

Short options are available for the input files and XTC output:

```bash
python3 center.py -s system.tpr -f trajectory.xtc \
  --selection 'protein and name CA' \
  --output-selection 'protein and name CA' \
  -o protein_ca_centered.xtc \
  --output-gro protein_ca_final.gro
```

If `--selection` is omitted, the script centers on `protein`. If `--output-selection` is omitted, it defaults to the centering selection. If `--output-gro` is provided, the script writes a GRO structure for the output selection from the final processed trajectory frame.

## Output notes

- The XTC contains only atoms matched by `--output-selection`.
- Pair the XTC with a structure containing the same output-selected atoms in the same order. When using `--output-gro`, use that GRO with the XTC.
- The XTC writer uses MDAnalysis precision setting `3`.
- The script does not modify the input files.

## References

- [MDAnalysis XTC reader and writer](https://docs.mdanalysis.org/stable/documentation_pages/coordinates/XTC.html)
- [MDAnalysis trajectory transformations](https://docs.mdanalysis.org/stable/documentation_pages/trajectory_transformations.html)
