"""Calculate successful divisions per actual cell update, grouped by state.

Run after a new simulation, for example:
    python calculate_birthrate.py data/log_eden_updates.tsv

World writes one update-log row for every cell whose action is processed. Its
state is measured immediately before that action. Deaths, skips, and blocked
birth attempts all contribute to the denominator; a successful division counts
once, regardless of the two resulting daughters. Final survivor snapshots are
not updates. This measures realized division frequency, not the configured
birth-attempt probability or a rate per hour.

The old lifetime log and subsampled ancestry logs lack the required update
history. Rerun the simulation with the updated World class to create it.
"""

#!/usr/bin/env python

import argparse
import csv
from pathlib import Path
import sys


DEFAULT_LOG = Path(__file__).resolve().parent / "data" / "log_eden_updates.tsv"


def calculate_birthrates(path: str | Path) -> dict[str, dict[str, int | float | None]]:
    """Return update counts, division counts, and division/update for each state.

    ``path`` must be a complete, unfiltered per-update TSV from one simulation.
    A state with no observed updates has a rate of None rather than zero.
    """
    counts = {state: {"cell_updates": 0, "divisions": 0} for state in ("edge", "center")}
    with open(path, newline="") as file:
        reader = csv.DictReader(file, delimiter="\t")
        required_columns = {"step", "id", "states", "action", "divided"}
        if not required_columns.issubset(reader.fieldnames or ()):
            raise ValueError(
                "Expected a per-update TSV with step, id, states, action, and divided columns. "
                "Rerun the simulation and use World.update_log_file, not the lifetime or subsampled log."
            )
        for row in reader:
            state = row["states"]
            action = row["action"]
            divided = row["divided"]
            if state not in counts:
                raise ValueError(f"Line {reader.line_num}: states must be 'edge' or 'center'.")
            if action not in ("birth", "death", "skip") or divided not in ("0", "1"):
                raise ValueError(f"Line {reader.line_num}: invalid action or divided flag.")
            if divided == "1" and action != "birth":
                raise ValueError(f"Line {reader.line_num}: only a birth action can divide.")
            counts[state]["cell_updates"] += 1
            counts[state]["divisions"] += int(divided)

    return {
        state: {
            **values,
            "birth_rate": values["divisions"] / values["cell_updates"]
            if values["cell_updates"] else None,
        }
        for state, values in counts.items()
    }


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description="Calculate successful divisions / actual cell updates for edge and center cells."
    )
    parser.add_argument(
        "update_log", nargs="?", type=Path, default=DEFAULT_LOG,
        help="complete *_updates.tsv from the updated simulator (default: %(default)s)",
    )
    args = parser.parse_args(argv)
    try:
        results = calculate_birthrates(args.update_log)
    except (OSError, ValueError) as error:
        parser.error(str(error))

    writer = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
    writer.writerow(("state", "cell_updates", "divisions", "birth_rate_per_update"))
    for state, values in results.items():
        rate = values["birth_rate"]
        writer.writerow((
            state, values["cell_updates"], values["divisions"],
            "NA" if rate is None else f"{rate:.8f}",
        ))


if __name__ == "__main__":
    main()
