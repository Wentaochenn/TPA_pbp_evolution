#!/usr/bin/env python3
"""Extract annotated Nichols PBP CDS translations, including initiation residues."""
import argparse
import hashlib
import json
from pathlib import Path
from Bio import SeqIO

GENES = {
    "mrcA": ("TPANIC_RS03500", 772432, 775087, -1, 884),
    "pbp1": ("TPANIC_RS02435", 535880, 537755, 1, 624),
    "pbp2": ("TPANIC_RS03760", 823826, 825785, 1, 652),
    "tp47": ("TPANIC_RS02830", 623857, 625162, -1, 434),
}


def extract(reference, genbank, output, audit=None):
    reference, genbank, output = map(Path, (reference, genbank, output))
    with genbank.open() as handle:
        records = list(SeqIO.parse(handle, "genbank"))
    if len(records) != 1 or records[0].id != "NC_021490.2":
        raise ValueError("Expected a single complete Nichols NC_021490.2 GenBank record")
    record = records[0]
    with reference.open() as handle:
        refs = list(SeqIO.parse(handle, "fasta"))
    if len(refs) != 1 or refs[0].id != record.id or str(refs[0].seq).upper() != str(record.seq).upper():
        raise ValueError("Reference FASTA and GenBank chromosome must agree exactly")
    sequences, rows = {}, []
    for gene, (tag, start, end, strand, length) in GENES.items():
        features = [f for f in record.features if f.type == "CDS" and
                    tag in f.qualifiers.get("locus_tag", [])]
        if len(features) != 1:
            raise ValueError("Expected exactly one annotated CDS for " + gene)
        feature = features[0]
        if (int(feature.location.start), int(feature.location.end), feature.location.strand) != (start, end, strand):
            raise ValueError("Reference CDS coordinates differ for " + gene)
        if feature.qualifiers.get("codon_start", ["1"]) != ["1"] or "transl_except" in feature.qualifiers:
            raise ValueError("Unsupported partial or exceptional CDS for " + gene)
        table = int(feature.qualifiers.get("transl_table", ["11"])[0])
        if table != 11:
            raise ValueError("Expected bacterial genetic code 11 for " + gene)
        cds = feature.extract(record.seq)
        protein = str(cds.translate(table=table, cds=True))
        annotated = feature.qualifiers.get("translation", [])
        if len(annotated) != 1 or protein != annotated[0] or len(protein) != length:
            raise ValueError("CDS translation does not match the annotated protein for " + gene)
        if not protein.startswith("M") or set(protein) - set("ACDEFGHIKLMNPQRSTVWY"):
            raise ValueError("Incomplete or unresolved reference protein for " + gene)
        sequences[gene] = protein
        rows.append(dict(gene=gene, locus_tag=tag, start_codon=str(cds[:3]),
                         first_residue=protein[0], length=length, translation_table=table,
                         matches_genbank_translation=True,
                         protein_sha256=hashlib.sha256(protein.encode()).hexdigest()))
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("".join(">" + g + "\n" + "\n".join(s[i:i+60] for i in range(0, len(s), 60)) + "\n"
                              for g, s in sequences.items()))
    if audit:
        result = dict(state="PASS", accession=record.id, proteins=rows,
                      reference_sha256=hashlib.sha256(reference.read_bytes()).hexdigest(),
                      genbank_sha256=hashlib.sha256(genbank.read_bytes()).hexdigest(),
                      protein_fasta_sha256=hashlib.sha256(output.read_bytes()).hexdigest())
        Path(audit).write_text(json.dumps(result, indent=2) + "\n")
    return sequences


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path, default=Path("refs/TPA_Nichols.fasta"))
    parser.add_argument("--genbank", type=Path, default=Path("refs/TPA_Nichols.gbff"))
    parser.add_argument("--output", type=Path, default=Path("data/nichols_pbp_proteins.fasta"))
    parser.add_argument("--audit", type=Path)
    args = parser.parse_args()
    proteins = extract(args.reference, args.genbank, args.output, args.audit)
    print("Wrote %d annotated Nichols PBP proteins to %s" % (len(proteins), args.output))


if __name__ == "__main__":
    main()
