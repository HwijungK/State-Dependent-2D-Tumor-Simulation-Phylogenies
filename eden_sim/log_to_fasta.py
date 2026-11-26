"""Convert an Eden simulation TSV to FASTA.

Usage:
    python eden_sim/log_to_fasta.py run.tsv --output run.fasta
"""

import argparse
from pathlib import Path

import pandas as pd


def convert_fasta(input_filename: str | Path, output_filename: str | Path) -> Path:
    """Write all logged cells to FASTA with state_id_death_time headers.

    Create output directories as needed, overwrite an existing output file,
    and return its path. Records are sorted by numeric cell ID.
    """
    df = pd.read_csv(
        input_filename,
        sep="\t",
        dtype={"id": "string", "death_date": "int64", "states": "string", "seq": "string"},
    )
    missing_columns = {"id", "death_date", "states", "seq"} - set(df.columns)
    if missing_columns:
        raise ValueError(f"Missing required columns: {', '.join(sorted(missing_columns))}.")

    df = df.sort_values(by="id", key=pd.to_numeric)

    output_filename = Path(output_filename)
    output_filename.parent.mkdir(parents=True, exist_ok=True)
    with output_filename.open("w") as fasta:
        for cell in df.itertuples(index=False):
            fasta.write(f">{cell.states}_{cell.id}_{cell.death_date}\n")
            fasta.write(f"{cell.seq}\n")
    return output_filename


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description="Convert all cells in an Eden simulation TSV to FASTA."
    )
    parser.add_argument("input", type=Path, help="simulation TSV containing id, death_date, states, and seq")
    parser.add_argument("-o", "--output", type=Path, required=True, help="output FASTA file")
    args = parser.parse_args(argv)
    try:
        output = convert_fasta(args.input, args.output)
    except (OSError, ValueError) as error:
        parser.error(str(error))
    print(f"Wrote {output}")


if __name__ == "__main__":
    main()
