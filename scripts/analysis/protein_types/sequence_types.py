from pathlib import Path
from collections import Counter, defaultdict
import csv, hashlib
BASES = frozenset("ACGT")
ALLOWED = frozenset("ACGTRYSWKMBDHVN?-")

def _write(path, columns, rows):
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)

def read_fasta(path):
    """Read a DNA alignment, retaining exact first-token FASTA identifiers."""
    sequences = {}
    current = None
    with Path(path).open(encoding="utf-8-sig") as handle:
        for number, raw in enumerate(handle, 1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                words = line[1:].split()
                if not words:
                    raise ValueError(f"Empty FASTA identifier at line {number}")
                current = words[0]
                if current in sequences:
                    raise ValueError(f"Duplicate FASTA identifier: {current}")
                sequences[current] = []
            else:
                if current is None:
                    raise ValueError(f"Sequence before FASTA header at line {number}")
                value = line.upper()
                illegal = set(value) - ALLOWED
                if illegal:
                    raise ValueError(f"Invalid DNA character(s) at line {number}: {sorted(illegal)}")
                sequences[current].append(value)
    sequences = {sample: "".join(parts) for sample, parts in sequences.items()}
    if len(sequences) < 2:
        raise ValueError("Regional alignment must contain at least two sequences")
    lengths = {len(seq) for seq in sequences.values()}
    if len(lengths) != 1 or not next(iter(lengths)):
        raise ValueError("FASTA sequences must have equal, nonzero aligned lengths")
    return sequences

def _group(metadata, sample, field):
    return str(metadata.get(sample, {}).get(field) or "UNKNOWN")

def _coordinates(settings, sequences):
    if settings.get("coordinate_system") != "reference_1based_inclusive_ungapped":
        raise ValueError("coordinate_system must be reference_1based_inclusive_ungapped")
    reference = settings.get("reference_sample")
    if reference not in sequences:
        raise ValueError("reference_sample must identify an entry in the regional FASTA")
    if set(sequences[reference]) - BASES:
        raise ValueError("Reference must contain only A/C/G/T; gapped or uncertain reference coordinates are unsupported")
    values = {}
    for field in ("region_start", "gene_start", "gene_end"):
        value = settings.get(field)
        if isinstance(value, bool) or not isinstance(value, int) or value < 1:
            raise ValueError(f"{field} must be a positive integer reference coordinate")
        values[field] = value
    end = values["region_start"] + len(sequences[reference]) - 1
    if not values["region_start"] <= values["gene_start"] <= values["gene_end"] <= end:
        raise ValueError("gene_start/gene_end must be an ordered interval inside the supplied region")
    # The full alignment is in forward reference coordinates. Strand is explicit
    # metadata; no amino-acid prediction or codon numbering is performed here.
    strand = settings.get("gene_strand")
    if strand not in ("+", "-"):
        raise ValueError("gene_strand must be explicitly + or -")
    return {**values, "region_end": end, "reference_sample": reference,
            "gene_strand": strand, "gene_slice": slice(values["gene_start"] - values["region_start"],
                                                       values["gene_end"] - values["region_start"] + 1)}

def _sequence_outputs(sequences, matrix, metadata, coords, outdir):
    gene_sequences = {sample: seq[coords["gene_slice"]] for sample, seq in sequences.items()}
    gene_length = len(next(iter(gene_sequences.values())))
    haplotypes = {}
    for seq in sorted(set(gene_sequences.values())):
        if set(seq) <= BASES:
            haplotypes[seq] = "H_" + hashlib.sha256(seq.encode("ascii")).hexdigest()[:16]
    rows = []
    group_counts = defaultdict(Counter)
    hap_counts = Counter()
    for sample in matrix:
        for field in ("lineage", "country", "source"):
            group_counts[(field, _group(metadata, sample, field))]["n_matrix_total"] += 1
    for sample, seq in sorted(sequences.items()):
        gene = gene_sequences[sample]
        n_gene_called = sum(base in BASES for base in gene)
        complete = n_gene_called == gene_length
        haplotype = haplotypes.get(gene, "UNASSIGNED_INCOMPLETE_SEQUENCE")
        lineage = _group(metadata, sample, "lineage")
        row = {"sample": sample, "lineage": lineage,
               "country": _group(metadata, sample, "country"),
               "source": _group(metadata, sample, "source"),
               "in_matrix": sample in matrix,
               "region_length": len(seq), "region_acgt_called": sum(base in BASES for base in seq),
               "gene_length": gene_length, "gene_acgt_called": n_gene_called,
               "gene_callability": n_gene_called / gene_length,
               "gene_complete_acgt": complete, "gene_haplotype": haplotype}
        rows.append(row)
        for field in ("lineage", "country", "source"):
            group_counts[(field, row[field])]["n_alignment"] += 1
            group_counts[(field, row[field])]["n_gene_complete"] += int(complete)
            group_counts[(field, row[field])]["n_incomplete"] += int(not complete)
        if complete:
            hap_counts[(haplotype, lineage)] += 1
    _write(outdir / "sequence_callability_and_haplotypes.tsv", list(rows[0]), rows)
    _write(outdir / "complete_case_selection_by_metadata.tsv",
           ["stratum_field", "stratum", "n_matrix_total", "n_alignment", "n_without_alignment",
            "n_gene_complete", "n_incomplete", "alignment_available_fraction", "complete_fraction",
            "complete_among_all_matrix_fraction"],
           [{"stratum_field": key[0], "stratum": key[1],
             **{field: counts[field] for field in ("n_matrix_total", "n_alignment", "n_gene_complete", "n_incomplete")},
             "n_without_alignment": counts["n_matrix_total"] - counts["n_alignment"],
             "alignment_available_fraction": counts["n_alignment"] / counts["n_matrix_total"],
             "complete_fraction": counts["n_gene_complete"] / counts["n_alignment"] if counts["n_alignment"] else None,
             "complete_among_all_matrix_fraction": counts["n_gene_complete"] / counts["n_matrix_total"]}
            for key, counts in sorted(group_counts.items())])
    _write(outdir / "gene_haplotype_lineage_counts.tsv", ["gene_haplotype", "lineage", "n_samples"],
           [{"gene_haplotype": key[0], "lineage": key[1], "n_samples": count}
            for key, count in sorted(hap_counts.items())])
    _write(outdir / "gene_haplotype_sequences.tsv", ["gene_haplotype", "length", "sequence_forward_reference_orientation"],
           [{"gene_haplotype": haplotype, "length": len(seq), "sequence_forward_reference_orientation": seq}
            for seq, haplotype in sorted(haplotypes.items(), key=lambda item: item[1])])
    alignment_samples = set(sequences)
    _write(outdir / "sample_scope_audit.tsv", ["sample", "scope"],
           [{"sample": sample, "scope": "matrix_without_regional_sequence"}
            for sample in sorted(set(matrix) - alignment_samples)] +
           [{"sample": sample, "scope": "regional_sequence_without_matrix_entry"}
            for sample in sorted(alignment_samples - set(matrix))])
    with (outdir / "gene_complete.fasta").open("w", encoding="utf-8") as handle:
        for sample, gene in sorted(gene_sequences.items()):
            if set(gene) <= BASES:
                handle.write(f">{sample}\n{gene}\n")
    with (outdir / "gene_haplotypes.fasta").open("w", encoding="utf-8") as handle:
        for seq, haplotype in sorted(haplotypes.items(), key=lambda item: item[1]):
            handle.write(f">{haplotype}\n{seq}\n")
    return gene_sequences, haplotypes, sum(row["gene_complete_acgt"] for row in rows)
