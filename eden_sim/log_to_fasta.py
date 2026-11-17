import pandas as pd

def convert_fasta(input_filename: str, output_filename: str):
    '''Convert log csv from simulation into a fasta file'''
    df = pd.read_csv(
        input_filename, sep = '\t',
        dtype = {"id": "string", "parent_id": "string", "birth_date": "int", "death_date": "int", "seq": "string"}
    )

    # Sort IDs in increasing order
    df = df.sort_values(by = "id", key=pd.to_numeric)

    # creates file if needed/overwrites it if alr exists
    with open(output_filename, "w") as fasta:
        # one FASTA entry for each logged cell
        for cell in df.itertuples(index=False):
            # build cell's dna sequence using its mutation positions

            # fasta file formatting
            fasta.write(f">cell_{cell.id}\n")
            fasta.write(f"{cell.seq}\n")