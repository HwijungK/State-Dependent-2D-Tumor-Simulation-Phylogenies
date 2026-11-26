#!/usr/bin/env python

"""Subsample final survivors and write their ancestry and end nodes as TSVs.

Usage:
    python bin/subsample.py --input run.tsv --output sampled.tsv --endnodes survivors.tsv --N 100 --seed 42
"""

import argparse
from numbers import Integral
from pathlib import Path

import pandas as pd


def subsample_log(
    path: str | Path,
    *,
    n: int,
    seed: int | None = None,
) -> pd.DataFrame:
    """Return sampled final survivors and every ancestor back to the origin.

    ``path`` points to a simulation TSV; ``n`` is a positive integer count.
    Sampling is uniform without replacement; ``seed`` makes it reproducible.
    Original columns, cell IDs, and row order are preserved.

    This uses the simulator's convention that final survivors have the maximum
    death_date, so the input must be a completed run with surviving cells.
    The returned row count includes ancestors and can exceed the sample size.
    """
    if isinstance(n, bool) or not isinstance(n, Integral) or n < 1:
        raise ValueError("n must be a positive integer.")

    df = pd.read_csv(
        Path(path),
        sep="\t",
        dtype={"id": "string", "parent_id": "string", "death_date": "int64"},
    )
    required_columns = {"id", "parent_id", "death_date"}
    missing_columns = required_columns - set(df.columns)
    if missing_columns:
        raise ValueError(f"Missing required columns: {', '.join(sorted(missing_columns))}.")
    if df.empty:
        raise ValueError("The simulation log is empty.")
    if df["id"].isna().any() or not df["id"].is_unique:
        raise ValueError("Each cell id must be unique and nonmissing.")
    if df["parent_id"].isna().sum() != 1:
        raise ValueError("Expected exactly one origin cell with a blank parent_id.")

    parents = df.set_index("id")["parent_id"].to_dict()
    missing_parents = set(df["parent_id"].dropna()) - set(parents)
    if missing_parents:
        raise ValueError(f"Missing parent IDs: {', '.join(sorted(missing_parents))}.")

    survivors = df.loc[df["death_date"].eq(df["death_date"].max()), "id"]
    sample_size = int(n)
    if sample_size > len(survivors):
        raise ValueError(
            f"Requested {sample_size} cells, but only {len(survivors)} final survivors are available."
        )
    sampled_ids = survivors.sample(n=sample_size, replace=False, random_state=seed)

    retained_ids = set()
    for cell_id in sampled_ids:
        lineage = set()
        while pd.notna(cell_id) and cell_id not in retained_ids:
            if cell_id in lineage:
                raise ValueError(f"Cycle detected in ancestry at cell {cell_id}.")
            lineage.add(cell_id)
            cell_id = parents[cell_id]
        # Only reuse a lineage after it has reached the origin or a verified path.
        retained_ids.update(lineage)

    return df.loc[df["id"].isin(retained_ids)].copy()


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description="Sample N final survivors, retaining their ancestry and exporting the end nodes separately."
    )
    parser.add_argument("-i", "--input", type=Path, required=True, help="complete simulation TSV")
    parser.add_argument("-o", "--output", type=Path, required=True, help="TSV of sampled survivors and ancestors")
    parser.add_argument("--endnodes", type=Path, required=True, help="TSV containing only the N sampled final survivors")
    parser.add_argument("-n", "--N", dest="n", type=int, required=True, help="positive integer number of survivors to sample")
    parser.add_argument("--seed", type=int, help="optional random seed for reproducible sampling")
    args = parser.parse_args(argv)
    try:
        paths = (args.input, args.output, args.endnodes)
        if len({path.resolve() for path in paths}) != len(paths):
            raise ValueError("input, output, and endnodes must be different files.")
        sampled = subsample_log(args.input, n=args.n, seed=args.seed)
        endnodes = sampled.loc[sampled["death_date"].eq(sampled["death_date"].max())]
        for path in (args.output, args.endnodes):
            path.parent.mkdir(parents=True, exist_ok=True)
        sampled.to_csv(args.output, sep="\t", index=False)
        endnodes.to_csv(args.endnodes, sep="\t", index=False)
    except (OSError, ValueError) as error:
        parser.error(str(error))
    print(f"Wrote {len(sampled)} survivors and ancestors to {args.output}")
    print(f"Wrote {len(endnodes)} final survivors to {args.endnodes}")


if __name__ == "__main__":
    main()
