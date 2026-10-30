#!/usr/bin/env python3
"""Submit one Snakemake job through SGE with a host free-memory threshold."""

import argparse
import json
import re
import subprocess
from pathlib import Path


def properties(path):
    with Path(path).open() as handle:
        for _ in range(20):
            line = handle.readline()
            if line.startswith("# properties = "):
                return json.loads(line[len("# properties = "):])
            if not line:
                break
    raise ValueError("Snakemake jobscript has no properties header: " + str(path))


def command(jobscript, project, qsub, queue):
    prop = properties(jobscript)
    threads = max(1, int(prop.get("threads", 1)))
    resources = prop.get("resources", {})
    total_mb = int(resources.get("mem_mb", 8000))
    scheduler_slots = max(threads, int(resources.get("sge_slots", threads)))
    host_free_mb = max(1024, total_mb)
    rule = re.sub(r"[^A-Za-z0-9_]", "_", prop.get("rule", "task"))[:32]
    logs = Path(project) / "logs" / "scheduler"
    logs.mkdir(parents=True, exist_ok=True)
    return [str(qsub), "-sync", "y", "-terse", "-cwd", "-V", "-S", "/bin/bash", "-q", queue,
            "-N", "tp_" + rule, "-pe", "smp", str(scheduler_slots),
            "-l", "mem_free=%dM" % host_free_mb,
            "-o", str(logs), "-e", str(logs), str(jobscript)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", required=True)
    parser.add_argument("--queue", default="all.q")
    parser.add_argument("--qsub", default="qsub")
    parser.add_argument("jobscript")
    args = parser.parse_args()
    cmd = command(args.jobscript, args.project, args.qsub, args.queue)
    completed = subprocess.run(cmd, capture_output=True, text=True)
    if completed.returncode:
        raise SystemExit("SGE submission/job failed: %s\n%s" % (" ".join(cmd), completed.stderr.strip()))
    print(completed.stdout.strip(), flush=True)


if __name__ == "__main__":
    main()
