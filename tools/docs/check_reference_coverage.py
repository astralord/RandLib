#!/usr/bin/env python3

import re
import sys
from pathlib import Path


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    registry = (root / "tests/support/cases.inc").read_text(encoding="utf-8")
    classes = set(re.findall(r"\b[A-Z][A-Za-z0-9_]*Rand\b", registry))

    distribution_pages = root / "docs/reference/distributions"
    reference = "\n".join(
        path.read_text(encoding="utf-8")
        for path in sorted(distribution_pages.glob("*.md"))
    )
    documented = set(re.findall(r"\b[A-Z][A-Za-z0-9_]*Rand\b", reference))

    missing = sorted(classes - documented)
    if missing:
        print("Distribution classes missing from reference headings:", file=sys.stderr)
        for name in missing:
            print(f"  {name}", file=sys.stderr)
        return 1

    print(f"Mathematical reference covers all {len(classes)} registered distribution names.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
