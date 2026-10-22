#!/usr/bin/env python3
"""Verify the portable staging bundle without fitting models or contacting GitHub."""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import re
import sys

ROOT = Path(__file__).resolve().parents[1]
MODELS = ("with_intercept", "without_intercept")
CATEGORIES = {"R", "simulations", "plots", "docs", "tests"}
MANIFEST = ROOT / "PACKAGE_MANIFEST.csv"
EXCLUDED = {"PACKAGE_MANIFEST.csv", "VALIDATION.json"}


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def inventory():
    return [p for p in sorted(ROOT.rglob("*"))
            if p.is_file() and p.relative_to(ROOT).as_posix() not in EXCLUDED]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-manifest", action="store_true",
                        help="Rebuild the manifest after intentional package edits.")
    args = parser.parse_args()
    errors, warnings = [], []
    directories = {}
    for model in MODELS:
        path = ROOT / model
        directories[model] = sorted(p.name for p in path.iterdir() if p.is_dir())
        if set(directories[model]) != CATEGORIES:
            errors.append(f"Asymmetric model categories: {model}: {directories[model]}")
        for required in ["README.md"]:
            if not (path / required).is_file():
                errors.append(f"Missing entry: {model}/{required}")

    files = inventory()
    for path in ROOT.rglob("*"):
        rel = path.relative_to(ROOT).as_posix()
        if path.is_symlink():
            errors.append(f"Nonportable symlink: {rel}")
        if any(part in {".git", "__pycache__", "cache", "runtime", ".Rproj.user"} for part in path.parts[len(ROOT.parts):]):
            errors.append(f"Runtime/internal directory: {rel}")
    for path in files:
        rel = path.relative_to(ROOT).as_posix()
        if path.name in {".DS_Store", ".Renviron", ".Rhistory", ".RData", ".httr-oauth"}:
            errors.append(f"Session/private file: {rel}")
        if path.suffix in {".so", ".o", ".dylib", ".dll", ".log", ".aux", ".pyc"}:
            errors.append(f"Generated binary/build file: {rel}")
        if path.stat().st_size >= 100 * 1024 * 1024:
            errors.append(f"File exceeds 100 MiB: {rel}")
        elif path.stat().st_size >= 50 * 1024 * 1024:
            warnings.append(f"Large file: {rel}")
        if path.suffix.lower() in {".r", ".rmd", ".py", ".sh", ".cpp"}:
            content = path.read_text(encoding="utf-8")
            if re.search(r"/(?:Users|home)/[^\s'\"]+", content):
                errors.append(f"Machine-specific home path in code: {rel}")

    for md in ROOT.rglob("*.md"):
        if not md.exists():
            continue
        for link in re.findall(r"\]\(([^)]+)\)", md.read_text()):
            if re.match(r"(?:https?://|mailto:|#)", link):
                continue
            name = link.split("#", 1)[0]
            if name and not (md.parent / name).exists():
                errors.append(f"Broken navigation link: {md.relative_to(ROOT)} -> {link}")

    rows = [{"path": p.relative_to(ROOT).as_posix(), "bytes": p.stat().st_size,
             "sha256": sha256(p)} for p in files]
    if args.write_manifest and not errors:
        with MANIFEST.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=["path", "bytes", "sha256"])
            writer.writeheader()
            writer.writerows(rows)
    if not MANIFEST.exists():
        errors.append("Missing PACKAGE_MANIFEST.csv; run --write-manifest after preparation.")
    else:
        with MANIFEST.open(newline="", encoding="utf-8") as stream:
            expected = {r["path"]: r for r in csv.DictReader(stream)}
        actual = {r["path"]: r for r in rows}
        if set(expected) != set(actual):
            errors.append("Manifest file set differs from package inventory.")
        for name in sorted(set(expected) & set(actual)):
            if expected[name]["sha256"] != actual[name]["sha256"] or int(expected[name]["bytes"]) != actual[name]["bytes"]:
                errors.append(f"Manifest mismatch: {name}")
    result = {"status": "passed" if not errors else "failed", "scope": "package layout, portability, documentation links and file integrity",
              "model_categories": directories, "file_count": len(rows),
              "total_bytes": sum(r["bytes"] for r in rows),
              "errors": errors, "warnings": warnings,
              "exclusions_from_manifest": sorted(EXCLUDED),
              "note": "Model tests are recorded separately; no Monte Carlo study is launched by this check."}
    (ROOT / "VALIDATION.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
    return 0 if not errors else 1


if __name__ == "__main__":
    sys.exit(main())
