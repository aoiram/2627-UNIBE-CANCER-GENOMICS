#!/usr/bin/env python3
from pathlib import Path
import subprocess
import time

import yaml


PROJECT_DIR = Path(__file__).resolve().parent

with (PROJECT_DIR / "config_main.yaml").open() as yamlfile:
    config = yaml.safe_load(yamlfile)

project_dir = Path(config["WORKDIR"])
log_dir = project_dir / "logs" / "barplots"
log_dir.mkdir(parents=True, exist_ok=True)

log_path = log_dir / "log_running_time.txt"
start_time = time.time()

command = [
    "snakemake",
    "--snakefile",
    str(PROJECT_DIR / "workflow" / "mut-signatures.rules"),
    "--profile",
    config["SLURM_PROFILE"],
    "--configfile",
    str(PROJECT_DIR / "config_main.yaml"),
    "--directory",
    str(project_dir),
]

print(f"Start workflow on project: {config['PROJECT']}")

if config["BARPLOTS"]:
    print("Start channels bar plot production!")
    with log_path.open("a") as log_file:
        log_file.write(f"\nProject name: {config['PROJECT']}\n")
        log_file.write(f"Start time: {time.ctime()}\n")
        result = subprocess.run(command, stdout=log_file, stderr=subprocess.STDOUT, check=False)
        elapsed = time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time))
        log_file.write(f"Elapsed time: {elapsed}\n")
        log_file.write(f"Exit status: {result.returncode}\n")

    if result.returncode:
        raise SystemExit(f"Snakemake failed; see {log_path}")

    print("6 and 96 channels bar plot production is over!")
else:
    print("No need for channels plots.")

print("Pipeline has finished!")