#!/usr/bin/env python3
"""Main pipeline: runs each enabled step as its own Snakemake workflow.

To add a step: write a new .rules file in workflow/, add a config flag to
config_main.yaml, and add one entry to STEPS below.

Usage:
    ./main-pipeline.py                     # run every step enabled in the config
    ./main-pipeline.py --steps sigprofiler # run only the named step(s)
    ./main-pipeline.py --dry-run           # snakemake -n for each step
"""
import argparse
import subprocess
import time
from pathlib import Path

import yaml

PROJECT_DIR = Path(__file__).resolve().parent
CONFIG_FILE = PROJECT_DIR / "config_main.yaml"
WORKFLOW_DIR = PROJECT_DIR / "workflow"

# One entry per step: name (used for --steps and the log folder), config flag
# that enables it, snakefile to run, and a label for messages.
STEPS = [
    {
        "name": "barplots",
        "flag": "BARPLOTS",
        "snakefile": "mut-signatures.rules",
        "label": "6 and 96 channels bar plots",
    },
    {
        "name": "sigprofiler",
        "flag": "SIGPROF",
        "snakefile": "sigprofiler.rules",
        "label": "SigProfiler",
    },
]


def build_command(snakefile, config, dry_run):
    command = [
        "snakemake",
        "--snakefile", str(WORKFLOW_DIR / snakefile),
        "--profile", config["SLURM_PROFILE"],
        "--configfile", str(CONFIG_FILE),
        "--directory", str(config["WORKDIR"]),
    ]
    if dry_run:
        command.append("--dry-run")
    return command


def run_step(step, config, dry_run):
    """Run one snakefile, log its output, and stop the pipeline if it fails."""
    log_dir = Path(config["WORKDIR"]) / "logs" / step["name"]
    log_dir.mkdir(parents=True, exist_ok=True)
    log_path = log_dir / "log_running_time.txt"

    print(f"Start {step['label']}!")
    start = time.time()
    with log_path.open("a") as log_file:
        log_file.write(f"\nProject name: {config['PROJECT']}\n")
        log_file.write(f"Snakefile: {step['snakefile']}\n")
        log_file.write(f"Start time: {time.ctime()}\n")
        log_file.flush()
        result = subprocess.run(
            build_command(step["snakefile"], config, dry_run),
            stdout=log_file,
            stderr=subprocess.STDOUT,
            check=False,
        )
        elapsed = time.strftime("%H:%M:%S", time.gmtime(time.time() - start))
        log_file.write(f"Elapsed time: {elapsed}\n")
        log_file.write(f"Exit status: {result.returncode}\n")

    if result.returncode:
        raise SystemExit(f"{step['label']} failed; see {log_path}")
    print(f"{step['label']} is over! ({elapsed})")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--steps", nargs="+", choices=[s["name"] for s in STEPS],
                        help="Run only these steps, ignoring the config flags")
    parser.add_argument("--dry-run", action="store_true", help="Pass --dry-run to snakemake")
    args = parser.parse_args()

    with CONFIG_FILE.open() as yamlfile:
        config = yaml.safe_load(yamlfile)

    # Check all snakefiles exist before starting anything
    for step in STEPS:
        if not (WORKFLOW_DIR / step["snakefile"]).is_file():
            raise SystemExit(f"Snakefile not found: {WORKFLOW_DIR / step['snakefile']}")

    print(f"Start workflow on project: {config['PROJECT']}")
    pipeline_start = time.time()

    for step in STEPS:
        enabled = step["name"] in args.steps if args.steps else config.get(step["flag"], False)
        if enabled:
            run_step(step, config, args.dry_run)
        else:
            print(f"Skipping {step['label']}.")

    elapsed = time.strftime("%H:%M:%S", time.gmtime(time.time() - pipeline_start))
    print(f"Pipeline has finished! (total {elapsed})")


if __name__ == "__main__":
    main()
