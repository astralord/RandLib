#!/usr/bin/env python3
"""Extract one tagged release section from CHANGELOG.md."""

from __future__ import annotations

import argparse
import re
from pathlib import Path


def extract(changelog: str, tag: str) -> str:
    version = tag.removeprefix("v")
    heading = re.compile(rf"^## \[{re.escape(version)}\](?:\s+-.*)?$", re.MULTILINE)
    match = heading.search(changelog)
    if match is None:
        raise ValueError(f"CHANGELOG.md has no [{version}] section")
    next_heading = re.search(r"^## \[", changelog[match.end() :], re.MULTILINE)
    end = match.end() + next_heading.start() if next_heading else len(changelog)
    body = changelog[match.end() : end].strip()
    if not body:
        raise ValueError(f"CHANGELOG.md section [{version}] is empty")
    return body + "\n"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("tag")
    parser.add_argument("--changelog", type=Path, default=Path("CHANGELOG.md"))
    arguments = parser.parse_args()
    print(extract(arguments.changelog.read_text(encoding="utf-8"), arguments.tag), end="")


if __name__ == "__main__":
    main()
