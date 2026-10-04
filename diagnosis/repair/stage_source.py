#!/usr/bin/env python3
"""Snapshot the current package and checks into a fresh FACE-HD scratch directory."""

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    output = args.output.resolve()
    if Path("/scratch.global/zhan9381/FACE-HD") not in output.parents or output.exists():
        parser.error("Choose a new output directory below the FACE-HD scratch root")
    names = subprocess.check_output(["git", "ls-files", "-co", "--exclude-standard", "-z"])
    paths = {Path(name.decode()) for name in names.split(b"\0") if name}
    roots = {"R", "src", "inst", "man", "tests", "scripts", "docs"}
    top_files = {"DESCRIPTION", "NAMESPACE", "LICENSE", "README.md", "main.R", "main.sh", "main.cmd", "realdata.R",
                 ".Rbuildignore", ".Rinstignore"}
    selected = sorted(path for path in paths if path.is_file() and
                      (path.parts[0] in roots or str(path) in top_files or
                       str(path).startswith("diagnosis/repair/")))
    output.mkdir(parents=True)
    source = output / "source"
    manifest = []
    for path in selected:
        destination = source / path
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, destination)
        manifest.append(dict(path=str(path), sha256=hashlib.sha256(destination.read_bytes()).hexdigest()))
    (output / "source_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print("Snapshot of {} files: {}".format(len(selected), source))


if __name__ == "__main__":
    main()
