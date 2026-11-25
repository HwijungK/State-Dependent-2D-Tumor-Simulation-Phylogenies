from pathlib import Path
import pandas as pd

DATA_DIR = Path(__file__).resolve().parent / "data"


def convert_fasta(input_filename: str | Path = DATA_DIR / "log_eden.tsv",
                  output_filename: str | Path = DATA_DIR / "log.fasta"):
    '''Convert the simulation log to FASTA with state_id_death_time headers.

    Input and output default to eden_sim/data.
    '''
    df = pd.read_csv(
        input_filename, sep = '\t',
        dtype = {"id": "string", "parent_id": "string", "birth_date": "int", "death_date": "int", "states": "string", "seq": "string"}
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
            fasta.write(f">{cell.states}_{cell.id}_{cell.death_date}\n")
            fasta.write(f"{cell.seq}\n")
