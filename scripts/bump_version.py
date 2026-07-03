#!/usr/bin/env python3
"""Bump the ProteinDJ version across all relevant files.

Usage:
    python bump_version.py 2.2.1

Updates:
  - nextflow.config          manifest.version (full)  and  params.container_version (major.minor, prefixed v)
  - bindsweeper/pyproject.toml  version (full)
  - apptainer/*.def          %labels Version (major.minor, no prefix)
"""

import re
import sys
from pathlib import Path


def parse_version(version_str: str) -> tuple[str, str, str]:
    """Parse 'MAJOR.MINOR.PATCH' and return (major, minor, patch)."""
    match = re.fullmatch(r"(\d+)\.(\d+)\.(\d+)", version_str)
    if not match:
        sys.exit(f"Error: version must be in MAJOR.MINOR.PATCH format, got '{version_str}'")
    return match.group(1), match.group(2), match.group(3)


def replace_in_file(path: Path, pattern: str, replacement: str) -> bool:
    """Replace first match of pattern in file. Returns True if a change was made."""
    text = path.read_text()
    new_text, count = re.subn(pattern, replacement, text, count=1, flags=re.MULTILINE)
    if count == 0:
        print(f"  WARNING: pattern not found in {path}")
        return False
    if new_text == text:
        print(f"  (no change) {path}")
        return False
    path.write_text(new_text)
    return True


def main() -> None:
    if len(sys.argv) != 2:
        sys.exit(f"Usage: {sys.argv[0]} MAJOR.MINOR.PATCH")

    new_full = sys.argv[1]
    major, minor, patch = parse_version(new_full)
    new_major_minor = f"{major}.{minor}"          # e.g. 2.2
    new_container_version = f"v{new_major_minor}" # e.g. v2.2

    root = Path(__file__).parent.parent


    print(f"Bumping ProteinDJ to {new_full} (container version {new_container_version})\n")

    # ── nextflow.config ────────────────────────────────────────────────────────
    nf_config = root / "nextflow.config"
    print(f"[nextflow.config]")
    changed = replace_in_file(
        nf_config,
        r"(container_version\s*=\s*')[^']*(')",
        rf"\g<1>{new_container_version}\g<2>",
    )
    if changed:
        print(f"  container_version → '{new_container_version}'")
    changed = replace_in_file(
        nf_config,
        r"(^\s*version\s*=\s*')[^']*(')",
        rf"\g<1>{new_full}\g<2>",
    )
    if changed:
        print(f"  manifest.version  → '{new_full}'")

    # ── bindsweeper/pyproject.toml ─────────────────────────────────────────────
    pyproject = root / "bindsweeper" / "pyproject.toml"
    print(f"\n[bindsweeper/pyproject.toml]")
    changed = replace_in_file(
        pyproject,
        r'(^version\s*=\s*")[^"]*(")',
        rf"\g<1>{new_full}\g<2>",
    )
    if changed:
        print(f"  version → \"{new_full}\"")

    # ── apptainer/*.def ────────────────────────────────────────────────────────
    def_files = sorted((root / "apptainer").glob("*.def"))
    print(f"\n[apptainer/*.def]  (label Version → {new_major_minor})")
    for def_file in def_files:
        changed = replace_in_file(
            def_file,
            r"(^[ \t]*Version\s+)\S+",
            rf"\g<1>{new_major_minor}",
        )
        if changed:
            print(f"  {def_file.name}")

    print(f"\nDone. All files updated to {new_full}.")


if __name__ == "__main__":
    main()
