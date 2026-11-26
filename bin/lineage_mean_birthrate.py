#!/usr/bin/env python

"""Calculate and export lineage mean birth rates for final simulation survivors.

Usage from the project directory:
    python lineage_mean_birthrate.py
    python lineage_mean_birthrate.py run.tsv --output rates.tsv

Preserves the notebook's definition: generation / (max(death_date) - 1).
The founder starts at generation 1. Final survivors are logged by end_sim()
with death_date equal to the final simulation tick plus one.
"""

import argparse
import csv
from pathlib import Path


DATA_DIR = Path(__file__).resolve().parent / "data"
DEFAULT_LOG = DATA_DIR / "log_eden.tsv"
DEFAULT_OUTPUT = DATA_DIR / "survivor_lineage_mean_birthrate.tsv"


def calculate_lineage_mean_birthrate(path: str | Path) -> dict[str, float]:
    """Return cell ID to generation/tick for survivors of a completed run."""
    with open(path, newline="") as file:
        reader = csv.DictReader(file, delimiter="\t")
        missing_columns = {"id", "generation", "death_date"} - set(reader.fieldnames or ())
        if missing_columns:
            raise ValueError(f"Missing required columns: {', '.join(sorted(missing_columns))}.")
        cells = [
            (row["id"], int(row["generation"]), int(row["death_date"]))
            for row in reader
        ]

    if not cells:
        raise ValueError("The simulation log is empty.")
    final_death_date = max(death_date for _, _, death_date in cells)
    tick = final_death_date - 1
    if tick <= 0:
        raise ValueError(
            "lineage_mean_birthrate = generation / tick requires at least one completed simulation tick."
        )

    return {
        cell_id: generation / tick
        for cell_id, generation, death_date in cells
        if death_date == final_death_date
    }


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description="Export final survivors' generation / simulation tick as a TSV."
    )
    parser.add_argument(
        "log_file", nargs="?", type=Path, default=DEFAULT_LOG,
        help="completed simulation lifetime TSV (default: %(default)s)",
    )
    parser.add_argument(
        "-o", "--output", type=Path, default=DEFAULT_OUTPUT,
        help="output TSV file (default: %(default)s)",
    )
    args = parser.parse_args(argv)
    try:
        rates = calculate_lineage_mean_birthrate(args.log_file)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open("w", newline="") as file:
            writer = csv.writer(file, delimiter="\t", lineterminator="\n")
            writer.writerow(("id", "lineage_mean_birthrate"))
            writer.writerows(rates.items())
    except (OSError, ValueError) as error:
        parser.error(str(error))
    print(f"Wrote {len(rates)} survivors to {args.output}")


if __name__ == "__main__":
    main()
