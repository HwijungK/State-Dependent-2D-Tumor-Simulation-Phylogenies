"""Fill the SDEVO BEAST template from a subsampled simulation FASTA.

Examples, run from the project directory:
    python eden_sim/create_beast_xml.py eden_sim/data
    python eden_sim/create_beast_xml.py path/to/sample.fasta --simulation-log path/to/log_eden.tsv

FASTA headers must be state_id_death_time, as written by log_to_fasta.py.
Sequences, type traits, alignment IDs/references, sampling proportions, and
BEAST output filenames are replaced. Sampling proportions are calculated
separately for center and edge: sampled / final survivors in the full log.
Other model parameters, priors, operators, and MCMC settings stay as in the
template. The template's contemporaneous sampling model requires all sampled
cells to have the same death_time.

Generated XML defaults to eden_sim/data/<FASTA stem>_SDEVO.xml. Existing output
files and the template are never overwritten; use --output for another name.
"""

import argparse
from collections import Counter
import csv
from dataclasses import dataclass
from decimal import Decimal
from pathlib import Path
import re
import xml.etree.ElementTree as ET


DATA_DIR = Path(__file__).resolve().parent / "data"
DEFAULT_TEMPLATE = DATA_DIR / "sdevor2" / "subsampled30_SDEVO_fixed.xml"
TYPES = ("center", "edge")  # BDMM-Prime TypeSet sorts type names alphabetically.
FASTA_SUFFIXES = {".fasta", ".fa", ".fas", ".fna"}


@dataclass(frozen=True)
class FastaRecord:
    taxon: str
    state: str
    cell_id: str
    death_time: Decimal
    sequence: str


def find_fasta(path: str | Path) -> Path:
    """Accept a FASTA file or a directory containing one subsampled FASTA."""
    path = Path(path)
    if path.is_dir():
        candidates = sorted(p for p in path.iterdir() if p.is_file() and p.suffix.lower() in FASTA_SUFFIXES)
        subsampled = [p for p in candidates if "subsampled" in p.stem.lower()]
        candidates = subsampled or candidates
        if len(candidates) != 1:
            raise ValueError(
                "Directory must contain one subsampled FASTA; pass the specific file if ambiguous."
            )
        path = candidates[0]
    if not path.is_file():
        raise FileNotFoundError(f"FASTA file not found: {path}")
    return path


def read_fasta(path: str | Path) -> list[FastaRecord]:
    """Read aligned, possibly multiline DNA sequences and simulation headers."""
    sequences = {}
    current_taxon = None
    with open(path) as file:
        for line_number, line in enumerate(file, start=1):
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                header = line[1:].split()
                current_taxon = header[0] if header else ""
                if not current_taxon or current_taxon in sequences:
                    raise ValueError(f"Line {line_number}: missing or duplicate FASTA header.")
                sequences[current_taxon] = []
            else:
                if current_taxon is None:
                    raise ValueError(f"Line {line_number}: sequence appears before its FASTA header.")
                sequences[current_taxon].append("".join(line.split()).upper())

    records = []
    cell_ids = set()
    for taxon, chunks in sequences.items():
        match = re.fullmatch(r"(center|edge)_([0-9]+)_([0-9]+(?:\.[0-9]+)?)", taxon)
        if not match:
            raise ValueError(f"Invalid FASTA header {taxon!r}; expected state_id_death_time.")
        state, cell_id, death_time = match.groups()
        cell_id = str(int(cell_id))
        if cell_id in cell_ids:
            raise ValueError(f"Cell {cell_id} appears more than once in the FASTA.")
        cell_ids.add(cell_id)
        sequence = "".join(chunks)
        if not sequence or set(sequence) - set("ACGTRYSWKMBDHVN?-"):
            raise ValueError(f"Taxon {taxon!r} needs a nonempty DNA sequence.")
        records.append(FastaRecord(taxon, state, cell_id, Decimal(death_time), sequence))
    if not records or len({len(record.sequence) for record in records}) != 1:
        raise ValueError("FASTA must contain nonempty sequences of equal length.")
    if len({record.death_time for record in records}) != 1:
        raise ValueError("This template requires contemporaneous samples with the same death_time.")
    return records


def calculate_sampling_proportions(
    records: list[FastaRecord], simulation_log: str | Path,
) -> dict[str, float]:
    """Count final survivors in the complete log and verify the sampled cells."""
    latest_time = None
    survivors = {}
    all_ids = set()
    children = Counter()
    roots = 0
    with open(simulation_log, newline="") as file:
        reader = csv.DictReader(file, delimiter="\t")
        required = {"id", "parent_id", "death_date", "states"}
        if not required.issubset(reader.fieldnames or ()):
            raise ValueError("The full simulation TSV needs id, parent_id, death_date, and states columns.")
        for row in reader:
            cell_id = row["id"]
            if not cell_id or cell_id in all_ids or row["states"] not in TYPES:
                raise ValueError(f"Line {reader.line_num}: invalid cell ID or state in simulation log.")
            all_ids.add(cell_id)
            if row["parent_id"]:
                children[row["parent_id"]] += 1
            else:
                roots += 1
            death_time = Decimal(row["death_date"])
            if not death_time.is_finite() or death_time < 0:
                raise ValueError(f"Line {reader.line_num}: invalid death_date.")
            if latest_time is None or death_time > latest_time:
                latest_time = death_time
                survivors = {}
            if death_time == latest_time:
                survivors[cell_id] = (row["states"], row.get("seq"))
    if roots != 1 or set(children) - all_ids or any(count != 2 for count in children.values()):
        raise ValueError("Use the complete two-daughter simulation log, not a subsampled ancestry log.")

    for record in records:
        if record.death_time != latest_time or record.cell_id not in survivors:
            raise ValueError(f"Sample {record.taxon!r} is not a final survivor in the simulation log.")
        state, sequence = survivors[record.cell_id]
        if record.state != state:
            raise ValueError(f"Sample {record.taxon!r} disagrees with its logged state.")
        if sequence is not None and record.sequence != sequence.upper():
            raise ValueError(f"Sample {record.taxon!r} disagrees with its logged sequence.")
    sampled = Counter(record.state for record in records)
    totals = Counter(state for state, _ in survivors.values())
    if any(totals[state] == 0 for state in TYPES):
        raise ValueError("Cannot infer sampling proportions: each modeled state needs final survivors.")
    return {state: sampled[state] / totals[state] for state in TYPES}


def _validate_references(root: ET.Element) -> None:
    ids = [element.attrib["id"] for element in root.iter() if "id" in element.attrib]
    if len(ids) != len(set(ids)):
        raise ValueError("Generated XML contains duplicate IDs.")
    known_ids = set(ids)
    for element in root.iter():
        for attribute, value in element.attrib.items():
            reference = value if attribute == "idref" else value[1:] if value.startswith("@") else None
            if reference is not None and reference not in known_ids:
                raise ValueError(f"Generated XML references undefined ID {reference!r}.")


def create_beast_xml(
    fasta_path: str | Path = DATA_DIR / "log_eden_subsampled.fasta",
    *,
    simulation_log: str | Path | None = None,
    template_path: str | Path = DEFAULT_TEMPLATE,
    output_path: str | Path | None = None,
) -> Path:
    """Fill the template and return the path to a newly created XML file.

    For ``*_subsampled.fasta``, the full log defaults to the neighboring TSV
    without ``_subsampled``. For other names it defaults to ``log_eden.tsv``
    beside the FASTA; callers can supply ``simulation_log`` explicitly.
    """
    fasta_path = find_fasta(fasta_path)
    records = read_fasta(fasta_path)
    if simulation_log is None:
        log_stem = (
            fasta_path.stem.removesuffix("_subsampled")
            if fasta_path.stem.endswith("_subsampled") else "log_eden"
        )
        simulation_log = fasta_path.with_name(f"{log_stem}.tsv")
    if not Path(simulation_log).is_file():
        raise FileNotFoundError(f"Full simulation log not found: {simulation_log}. Supply --simulation-log.")
    proportions = calculate_sampling_proportions(records, simulation_log)
    template_path = Path(template_path)
    output_path = Path(output_path) if output_path is not None else DATA_DIR / f"{fasta_path.stem}_SDEVO.xml"
    if output_path.resolve() == template_path.resolve():
        raise ValueError("Output must be a new file, not the template.")
    if output_path.exists():
        raise FileExistsError(f"Output already exists: {output_path}. Choose another --output filename.")

    parser = ET.XMLParser(target=ET.TreeBuilder(insert_comments=True))
    tree = ET.parse(template_path, parser=parser)
    root = tree.getroot()
    alignments = root.findall("data")
    if len(alignments) != 1 or not alignments[0].get("id"):
        raise ValueError("Template must contain one named data alignment.")
    alignment = alignments[0]
    old_id = alignment.attrib["id"]
    alignment_id = re.sub(r"[^A-Za-z0-9_.-]", "_", fasta_path.stem)
    old_sequences = alignment.findall("sequence")
    if not old_sequences:
        raise ValueError("Template alignment has no sequence to use as a prototype.")
    sequence_attributes = dict(old_sequences[0].attrib)
    for sequence in old_sequences:
        alignment.remove(sequence)
    for record in records:
        attributes = {
            **sequence_attributes, "id": f"seq_{record.taxon}",
            "taxon": record.taxon, "value": record.sequence,
        }
        ET.SubElement(alignment, "sequence", attributes)

    type_traits = root.findall(".//typeTraitSet[@traitname='type']")
    if len(type_traits) != 1:
        raise ValueError("Template must contain one type trait mapping.")
    type_traits[0].set("value", ",".join(f"{record.taxon}={record.state}" for record in records))
    for type_set in root.findall(".//typeSet[@spec]"):
        # Keep both model states even when the sample contains just one state.
        type_set.set("value", ",".join(TYPES))
    rho_values = root.findall(".//rhoSampling/parameter[@name='values']")
    if len(rho_values) != 1:
        raise ValueError("Template must contain one rhoSampling values parameter.")
    rho_values[0].set("dimension", str(len(TYPES)))
    rho_values[0].text = " ".join(str(proportions[state]) for state in TYPES)

    for element in root.iter():
        for attribute, value in list(element.attrib.items()):
            element.set(attribute, value.replace(old_id, alignment_id))
        filename = element.get("fileName")
        if filename:
            for prefix in ("$(filebase).$(tree)", "$(filebase)", "$(tree)", alignment_id):
                if filename.startswith(prefix):
                    element.set("fileName", output_path.stem + filename[len(prefix):])
                    break

    _validate_references(root)
    ET.indent(tree, space="    ")
    xml = ET.tostring(root, encoding="utf-8", xml_declaration=True) + b"\n"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("xb") as file:
        file.write(xml)
    return output_path


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description="Create a BEAST SDEVO XML from subsampled simulation sequences."
    )
    parser.add_argument(
        "fasta", nargs="?", type=Path, default=DATA_DIR / "log_eden_subsampled.fasta",
        help="FASTA file or directory containing one subsampled FASTA",
    )
    parser.add_argument(
        "--simulation-log", type=Path,
        help="full simulation TSV used to calculate per-state sampling proportions",
    )
    parser.add_argument("--template", type=Path, default=DEFAULT_TEMPLATE, help="SDEVO XML template")
    parser.add_argument(
        "-o", "--output", type=Path,
        help="new XML filename (default: eden_sim/data/<FASTA stem>_SDEVO.xml)",
    )
    args = parser.parse_args(argv)
    try:
        output = create_beast_xml(
            args.fasta, simulation_log=args.simulation_log,
            template_path=args.template, output_path=args.output,
        )
    except (OSError, ValueError, ArithmeticError, ET.ParseError) as error:
        parser.error(str(error))
    print(f"Wrote {output}")


if __name__ == "__main__":
    main()
