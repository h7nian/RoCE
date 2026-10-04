#!/usr/bin/env python3
"""Build a hashed source manifest and manuscript claim index, without execution.

R functions are indexed separately using the R parser. Python functions use
the Python AST. C++ definitions are lexical candidates, explicitly not a full
C++ AST or a semantic review. Existing generated/result directories are excluded.
"""

import argparse
import ast
import csv
import hashlib
from pathlib import Path
import re
import subprocess


def write_csv(path, rows, fields):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    output = args.output.resolve()
    if Path("/scratch.global") not in output.parents or output.exists():
        parser.error("Choose a new output directory under /scratch.global")
    candidates = subprocess.check_output(["rg", "--files", "--hidden", "-g", "!.git/**"],
                                         universal_newlines=True).splitlines()
    suffixes = {".R", ".r", ".cpp", ".h", ".hpp", ".py", ".sh", ".cmd"}
    configuration_names = {"DESCRIPTION", "NAMESPACE", "Makevars", "Makefile", ".Rprofile"}
    excluded_parts = {"results", "out", "logs", "overleaf", ".git", ".Rproj.user"}
    paths = sorted(Path(p) for p in candidates if
                   (Path(p).suffix in suffixes or Path(p).name in configuration_names) and
                   not excluded_parts.intersection(Path(p).parts))
    # Rcpp registration wrappers are git-ignored but participate in the package
    # interface. Include the existing files without regenerating them.
    generated_interfaces = [Path("R/RcppExports.R"), Path("src/RcppExports.cpp")]
    paths = sorted(set(paths).union(path for path in generated_interfaces if path.is_file()))
    manifest, functions, parse_checks, claims = [], [], [], []
    for path in paths:
        data = path.read_bytes()
        source = data.decode("utf-8")
        manifest.append(dict(path=str(path), sha256=hashlib.sha256(data).hexdigest(),
                             bytes=len(data), lines=len(source.splitlines()),
                             semantic_review="pending"))
        if path.suffix == ".py":
            try:
                tree = ast.parse(source)
            except SyntaxError as error:
                parse_checks.append(dict(path=str(path), status="syntax_error", detail=str(error)))
                continue
            parse_checks.append(dict(path=str(path), status="parsed", detail="Python AST only"))
            for node in ast.walk(tree):
                if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.Lambda)):
                    functions.append(dict(path=str(path), line=node.lineno,
                                          name=getattr(node, "name", "<lambda>"),
                                          mechanism="python_ast", semantic_review="pending"))
        elif path.suffix in {".cpp", ".h", ".hpp"}:
            # Mask comments and strings while preserving all newline positions.
            masked = re.sub(r'//[^\n]*|/\*[\s\S]*?\*/|"(?:\\.|[^"\\])*"',
                            lambda match: re.sub(r"[^\n]", " ", match.group()), source)
            pattern = r"(?:[\w:<>,*&]+\s+)+([\w:]+)\s*\([^;{}]*?\)\s*(?:const\s*)?\{"
            for match in re.finditer(pattern, masked):
                if match.group(1) in {"if", "for", "while", "switch", "catch"}:
                    continue
                functions.append(dict(path=str(path), line=source.count("\n", 0, match.start()) + 1,
                                      name=match.group(1), mechanism="cpp_lexical_candidate",
                                      semantic_review="pending"))
        elif path.suffix in {".sh", ".cmd"}:
            for match in re.finditer(r"(?m)^\s*(?:function\s+)?(\w+)\s*\(\)\s*\{", source):
                functions.append(dict(path=str(path), line=source.count("\n", 0, match.start()) + 1,
                                      name=match.group(1), mechanism="shell_lexical_candidate",
                                      semantic_review="pending"))
    for path in [Path("docs/main.tex"), Path("docs/supplemental.tex")]:
        source = path.read_text()
        pattern = r"\\begin\{(assumption|lemma|theorem|proposition)\}(?:\[([^\]]*)\])?"
        for match in re.finditer(pattern, source):
            end = source.find("\\end{" + match.group(1) + "}", match.end())
            block = source[match.end():end]
            label = re.search(r"\\label\{([^}]+)\}", block)
            claims.append(dict(path=str(path), line=source.count("\n", 0, match.start()) + 1,
                               kind=match.group(1), title=match.group(2) or "",
                               label=label.group(1) if label else "",
                               review_reference="theory_report.md and claim_review.csv"))
    output.mkdir(parents=True)
    write_csv(output / "source_manifest.csv", manifest,
              ["path", "sha256", "bytes", "lines", "semantic_review"])
    write_csv(output / "non_r_function_index.csv", functions,
              ["path", "line", "name", "mechanism", "semantic_review"])
    write_csv(output / "python_parse_checks.csv", parse_checks, ["path", "status", "detail"])
    write_csv(output / "manuscript_claim_index.csv", claims,
              ["path", "line", "kind", "title", "label", "review_reference"])
    print("Indexed {} source files, {} non-R function entries, {} manuscript claims.".format(
        len(manifest), len(functions), len(claims)))


if __name__ == "__main__":
    main()
