#!/usr/bin/env python3
"""Create an isolated project root for the original tp_wgs analysis scripts."""

import argparse
import csv
import json
import shutil
from pathlib import Path


def link(source, target):
    if not source.exists():
        raise FileNotFoundError(str(source))
    if target.is_symlink():
        if target.resolve() != source.resolve():
            raise ValueError("Existing link points elsewhere: " + str(target))
    elif target.exists():
        raise ValueError("Refusing to overwrite existing input: " + str(target))
    else:
        target.parent.mkdir(parents=True, exist_ok=True)
        target.symlink_to(source.resolve())


def copy_input(source, target):
    if not source.is_file():
        raise FileNotFoundError(str(source))
    if target.exists():
        if target.read_bytes() != source.read_bytes():
            raise ValueError("Existing input differs: " + str(target))
    else:
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-project", type=Path, required=True)
    parser.add_argument("--run-root", type=Path, required=True)
    args = parser.parse_args()
    source = args.source_project.resolve()
    run = args.run_root.resolve()
    package = Path(__file__).resolve().parent.parent
    if run == source or run == package or run in (source / "raw_fastq", source / "consensus", source / "bam"):
        parser.error("Choose a distinct run directory outside input directories")
    for name in ("raw_fastq", "refs", "consensus"):
        if not (source / name).is_dir():
            raise FileNotFoundError(str(source / name))
    for name in ("TPA_master_metadata.xlsx", "sample_dates.tsv"):
        if not (source / name).is_file():
            raise FileNotFoundError(str(source / name))
    excludes = set((package / "remove_samples.txt").read_text().split())
    sample_lists = {}
    for group in ("ss14", "nichols"):
        original = (package / "list" / ("tpa." + group + ".sampleList.txt")).read_text().split()
        if len(set(original)) != len(original):
            raise ValueError("Duplicate IDs in " + group + " sample list")
        selected = [sample for sample in original if sample not in excludes]
        for sample in selected:
            for end in ("1", "2"):
                fq = source / "raw_fastq" / group / (sample + "_" + end + ".fastq.gz")
                if not fq.is_file():
                    raise FileNotFoundError(str(fq))
        sample_lists[group] = selected
    if set(sample_lists["ss14"]) & set(sample_lists["nichols"]):
        raise ValueError("Sample IDs overlap between reference groups")
    with (package / "config" / "complete_genome_manifest.tsv").open() as handle:
        assemblies = list(csv.DictReader(handle, delimiter="\t"))
    if len(assemblies) != 5 or len({row["sample_id"] for row in assemblies}) != 5:
        raise ValueError("Expected five distinct assembly genomes")
    for row in assemblies:
        path = source / row["fasta"]
        if not path.is_file():
            raise FileNotFoundError(str(path))
    run.mkdir(parents=True, exist_ok=True)
    if not run.is_dir():
        raise ValueError("Could not create run directory")
    for name in ("script", "SC", "envs"):
        link(package / name, run / name)
    link(source / "raw_fastq", run / "raw_fastq")
    for name in ("master.sh", "config.yaml", "remove_samples.txt"):
        copy_input(package / name, run / name)
    for name in ("TPA_master_metadata.xlsx", "sample_dates.tsv"):
        copy_input(source / name, run / name)
    for entry in (package / "config").iterdir():
        if entry.is_file():
            copy_input(entry, run / "config" / entry.name)
    for group, samples in sample_lists.items():
        dest = run / "list" / ("tpa." + group + ".sampleList.txt")
        dest.parent.mkdir(parents=True, exist_ok=True)
        content = "".join(sample + "\n" for sample in samples)
        if dest.exists() and dest.read_text() != content:
            raise ValueError("Filtered sample list differs: " + str(dest))
        dest.write_text(content)
    if not excludes.issubset(set((package / "list/tpa.ss14.sampleList.txt").read_text().split()) | set((package / "list/tpa.nichols.sampleList.txt").read_text().split())):
        raise ValueError("Excluded ID not found in original sample lists")
    for entry in (source / "refs").iterdir():
        if entry.is_file():
            copy_input(entry, run / "refs" / entry.name)
    for row in assemblies:
        link(source / row["fasta"], run / row["fasta"])
    pdb_dir = source / "data" / "alphafold_pdb"
    if pdb_dir.is_dir():
        dest = run / "data" / "alphafold_pdb"
        dest.mkdir(parents=True, exist_ok=True)
        for pdb in pdb_dir.glob("*.pdb"):
            link(pdb, dest / pdb.name)
    count = sum(map(len, sample_lists.values())) + len(assemblies)
    if count != 1861:
        raise ValueError("Expected 1861 selected genomes; found " + str(count))
    manifest = {"source_project": str(source), "run_root": str(run), "read_samples": {k: len(v) for k, v in sample_lists.items()}, "excluded_samples": sorted(excludes), "assembly_samples": [r["sample_id"] for r in assemblies], "total_selected": count}
    (run / "run_inputs.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print("Prepared %d selected genomes in %s" % (count, run))


if __name__ == "__main__":
    main()
