from pathlib import Path
import pandas as pd

DATA_DIR = Path(__file__).resolve().parent / "data"


def convert_fasta(input_filename: str | Path = DATA_DIR / "log_eden.tsv",
                  output_filename: str | Path = DATA_DIR / "log.fasta"):
    '''Convert the simulation log into a FASTA file, using eden_sim/data by default.'''
    df = pd.read_csv(
        input_filename, sep = '\t',
        dtype = {"id": "string", "parent_id": "string", "birth_date": "int", "death_date": "int", "seq": "string"}
    )

    # Sort IDs in increasing order
    df = df.sort_values(by = "id", key=pd.to_numeric)

    # creates the output directory and file if needed/overwrites it if alr exists
    output_filename = Path(output_filename)
    output_filename.parent.mkdir(parents = True, exist_ok = True)
    with open(output_filename, "w") as fasta:
        # one FASTA entry for each logged cell
        for cell in df.itertuples(index=False):
            # build cell's dna sequence using its mutation positions

            # fasta file formatting
            fasta.write(f">cell_{cell.id}\n")
            fasta.write(f"{cell.seq}\n")
