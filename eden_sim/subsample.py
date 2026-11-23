"""Subsample final simulation survivors while retaining their full ancestry."""

from math import ceil
from numbers import Integral, Real
from pathlib import Path

import pandas as pd


def subsample_log(
    path: str | Path,
    *,
    proportion: float | None = None,
    n: int | None = None,
    seed: int | None = None,
) -> pd.DataFrame:
    """Return sampled final survivors and every ancestor back to the origin.

    ``path`` points to a simulation TSV. Supply exactly one of ``proportion``
    (in (0, 1], rounded up to a whole cell) or ``n`` (a positive integer).
    Sampling is uniform without replacement; ``seed`` makes it reproducible.
    Original columns, cell IDs, and row order are preserved.

    This uses the simulator's convention that final survivors have the maximum
    death_date, so the input must be a completed run with surviving cells.
    The returned row count includes ancestors and can exceed the sample size.
    """
    if (proportion is None) == (n is None):
        raise ValueError("Supply exactly one of proportion or n.")

    if proportion is not None:
        if isinstance(proportion, bool) or not isinstance(proportion, Real) or not 0 < proportion <= 1:
            raise ValueError("proportion must be a number greater than 0 and at most 1.")
    elif isinstance(n, bool) or not isinstance(n, Integral) or n < 1:
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
    sample_size = ceil(proportion * len(survivors)) if proportion is not None else int(n)
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
