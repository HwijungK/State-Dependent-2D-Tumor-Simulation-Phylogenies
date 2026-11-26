"""Convert SDevo posterior birth/death rates to per-tick division probabilities.

For each MCMC sample and state, calculate
    p = birth / (birth + death) * (1 - exp(-(birth + death))).
This is the probability of division before death and before the tick ends,
assuming constant rates and no state changes within that tick.

The DataFrame retains the original columns and adds birth_probability_center,
birth_probability_edge, and center_minus_edge. The reported difference is the
mean of center_minus_edge across retained samples. All rows are used by
default; --burnin 0.1 discards the first 10 percent of logged samples.

Examples:
    python eden_sim/calculate-predicted-birthrate.py run.log
    python eden_sim/calculate-predicted-birthrate.py run.log --burnin 0.1 --output predictions.tsv

One model-time unit is assumed to equal one simulation tick for this study.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def calculate_predicted_birthrates(
    log_file: str | Path,
    *,
    burnin_proportion: float = 0.0,
) -> tuple[pd.DataFrame, float]:
    """Return the enhanced posterior DataFrame and mean(center - edge).

    Uses state-specific deathRateSVCanonical.center/edge when available,
    otherwise the shared deathRateSVCanonical column. Each posterior row is
    transformed separately, preserving correlations between birth and death.
    """
    if not np.isfinite(burnin_proportion) or not 0 <= burnin_proportion < 1:
        raise ValueError("burnin_proportion must be between 0 (inclusive) and 1 (exclusive).")

    df = pd.read_csv(log_file, sep="\t", comment="#")
    if df.empty:
        raise ValueError("The log contains no MCMC samples.")
    df = df.iloc[int(len(df) * burnin_proportion):].copy()

    for state in ("center", "edge"):
        birth_column = f"birthRateSVCanonical.{state}"
        death_column = "deathRateSVCanonical"

        missing = {birth_column, death_column} - set(df.columns)
        if missing:
            raise ValueError(f"Missing required rate columns: {', '.join(sorted(missing))}.")

        rates = {}
        for column in (birth_column, death_column):
            values = pd.to_numeric(df[column], errors="raise").to_numpy(dtype=float)
            if not np.isfinite(values).all() or (values < 0).any():
                raise ValueError(f"{column} must contain finite, nonnegative rates.")
            rates[column] = values
        birth = rates[birth_column]
        death = rates[death_column]
        total = birth + death

        # The zero-birth, zero-death case has probability zero. expm1 also
        # preserves accuracy when the total rate is close to zero.
        birth_share = np.divide(birth, total, out=np.zeros_like(birth), where=total > 0)
        df[f"birth_probability_{state}"] = birth_share * -np.expm1(-total)

    df["center_minus_edge"] = df["birth_probability_center"] - df["birth_probability_edge"]
    return df, float(df["center_minus_edge"].mean())


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description="Calculate per-tick birth probabilities and mean(center - edge) from a SDevo log."
    )
    parser.add_argument("log_file", type=Path, help="SDevo .log file")
    parser.add_argument(
        "--burnin", type=float, default=0.0,
        help="proportion of rows to discard as burn-in, from 0 to less than 1 (default: 0)",
    )
    parser.add_argument(
        "-o", "--output", type=Path,
        help="optionally save the enhanced DataFrame to a new TSV file",
    )
    args = parser.parse_args(argv)
    try:
        df, difference = calculate_predicted_birthrates(
            args.log_file, burnin_proportion=args.burnin,
        )
        if args.output is not None:
            args.output.parent.mkdir(parents=True, exist_ok=True)
            with args.output.open("x") as file:
                df.to_csv(file, sep="\t", index=False)
    except (OSError, ValueError) as error:
        parser.error(str(error))

    print(f"mean_probability_center\t{df['birth_probability_center'].mean():.12g}")
    print(f"mean_probability_edge\t{df['birth_probability_edge'].mean():.12g}")
    print(f"mean_center_minus_edge\t{difference:.12g}")
    if args.output is not None:
        print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
