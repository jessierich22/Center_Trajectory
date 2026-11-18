#!/usr/bin/env python3

import argparse
import warnings

warnings.filterwarnings(
    "ignore",
    message=r"The Bio\.Application modules and modules relying on it have been deprecated\.",
    module=r"Bio\.Application",
)

import MDAnalysis as mda
from MDAnalysis.exceptions import NoDataError

from MDAnalysis.transformations import (
    center_in_box,
    unwrap,
    wrap,
)


def main():

    parser = argparse.ArgumentParser(
        description=(
            "Reconstruct molecules across PBC, center an atom selection in "
            "the simulation box, and write a selected XTC trajectory."
        )
    )

    parser.add_argument(
        "-s",
        "--structure",
        required=True,
        help="Full-system topology with bond information (e.g. .tpr)",
    )

    parser.add_argument(
        "-f",
        "--trajectory",
        required=True,
        help="Full-system trajectory file (.xtc)",
    )

    parser.add_argument(
        "-o",
        "--output",
        default="centered.xtc",
        help="Output centered trajectory (.xtc)",
    )

    parser.add_argument(
        "--selection",
        default="protein",
        help=(
            "MDAnalysis atom selection used to center the system "
            "(default: 'protein')"
        ),
    )

    parser.add_argument(
        "--output-selection",
        default=None,
        help=(
            "MDAnalysis atom selection to write to the output trajectory "
            "(default: same as --selection)"
        ),
    )

    parser.add_argument(
        "--output-gro",
        default=None,
        help=(
            "Optional GRO structure for the output selection in the final "
            "trajectory frame"
        ),
    )

    args = parser.parse_args()

    print("=" * 70)
    print(" PROTEIN TRAJECTORY CENTERING")
    print("=" * 70)

    print(f"Structure        : {args.structure}")
    print(f"Trajectory       : {args.trajectory}")
    print(f"Output           : {args.output}")
    print(f"Center selection : {args.selection}")
    print(
        "Output selection : "
        f"{args.selection if args.output_selection is None else args.output_selection}"
    )
    if args.output_gro:
        print(f"Final GRO        : {args.output_gro}")
    print()

    # ============================================================
    # Load full system
    # ============================================================

    print("Loading trajectory...")

    universe = mda.Universe(
        args.structure,
        args.trajectory,
    )

    try:
        bonds = universe.bonds
    except NoDataError as exc:
        raise ValueError(
            "Molecule reconstruction requires bond information. "
            "Use a topology-bearing structure such as a GROMACS .tpr file; "
            "a .gro file alone does not contain bonds."
        ) from exc

    if len(bonds) == 0:
        raise ValueError(
            "Molecule reconstruction requires bond information, but no bonds "
            "were found in the input topology."
        )

    centering_atoms = universe.select_atoms(args.selection)

    if centering_atoms.n_atoms == 0:
        raise ValueError(
            f"The centering selection '{args.selection}' did not match any atoms."
        )

    output_selection = (
        args.selection
        if args.output_selection is None
        else args.output_selection
    )
    output_atoms = universe.select_atoms(output_selection)

    if output_atoms.n_atoms == 0:
        raise ValueError(
            f"The output selection '{output_selection}' did not match any atoms."
        )

    print(f"Total atoms      : {universe.atoms.n_atoms:,}")
    print(f"Center atoms     : {centering_atoms.n_atoms:,}")
    print(f"Output atoms     : {output_atoms.n_atoms:,}")
    print(f"Frames           : {len(universe.trajectory):,}")
    print()

    # ============================================================
    # PBC correction and centering
    # ============================================================

    print("Applying PBC correction and centering...")

    transformations = [
        unwrap(universe.atoms),
        center_in_box(
            centering_atoms,
            center="geometry",
            wrap=False,
        ),
        wrap(universe.atoms, compound="fragments"),
    ]

    universe.trajectory.add_transformations(*transformations)

    # ============================================================
    # Write selected atoms
    # ============================================================

    print("Writing centered selected-atoms trajectory...")

    with mda.Writer(
        args.output,
        output_atoms.n_atoms,
        precision=3,
    ) as writer:

        for timestep in universe.trajectory:
            writer.write(output_atoms)

    if args.output_gro:
        with mda.Writer(args.output_gro, output_atoms.n_atoms) as writer:
            writer.write(output_atoms)

    # ============================================================
    # Finished
    # ============================================================

    print()
    print("=" * 70)
    print(" Finished")
    print("=" * 70)
    print(f"Output        : {args.output}")
    print(f"Atoms written : {output_atoms.n_atoms:,}")
    if args.output_gro:
        print(f"Final GRO     : {args.output_gro}")
    print("XTC precision : 0.001 nm")
    print("=" * 70)


if __name__ == "__main__":
    main()
