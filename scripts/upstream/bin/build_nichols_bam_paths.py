#!/usr/bin/env python3
"""Write the selected Nichols-coordinate BAM paths for p2_02."""

import argparse
import csv
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root = args.project.resolve()
    paths = []
    for group, folder, suffix in (
        ("ss14", "ss14_align2_nichols", ".ss14_a2n.final.bam"),
        ("nichols", "nichols", ".nichols.final.bam"),
    ):
        samples = (root / "list" / ("tpa." + group + ".sampleList.txt")).read_text().split()
        paths.extend(root / "bam" / folder / "postfilter" / (s + suffix) for s in samples)
    with (root / "config" / "complete_genome_manifest.tsv").open() as handle:
        for record in csv.DictReader(handle, delimiter="\t"):
            paths.append(root / "results" / "04_complete_genome_variant_calls" / "bam" / (record["sample_id"] + "_vs_nichols.bam"))
    if len(paths) != 1861 or len(set(paths)) != 1861:
        raise ValueError("Expected 1861 distinct BAMs; found %d" % len(paths))
    missing = [p for p in paths if not p.is_file()]
    if missing:
        raise FileNotFoundError("Missing %d BAMs; first: %s" % (len(missing), missing[0]))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("".join(str(p) + "\n" for p in paths))
    print("Selected 1861 Nichols-coordinate BAMs: " + str(args.output))


if __name__ == "__main__":
    main()
