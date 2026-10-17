#!/usr/bin/env python3

import re
import sys
from pathlib import Path


def fail(errors: list[str], message: str) -> None:
    errors.append(message)


def main() -> int:
    if len(sys.argv) != 2:
        print(f"usage: {Path(sys.argv[0]).name} HTML_DIRECTORY", file=sys.stderr)
        return 2

    html_dir = Path(sys.argv[1])
    root = Path(__file__).resolve().parents[2]
    source_dir = root / "docs/reference/distributions"
    errors: list[str] = []
    labels: set[str] = set()
    card_labels: set[str] = set()

    for source in sorted(source_dir.glob("*.md")):
        text = source.read_text(encoding="utf-8")
        title = re.fullmatch(r"# (.+) \{#([a-z0-9_]+)\}", text.splitlines()[0])
        if title is None:
            fail(errors, f"{source.relative_to(root)}: missing a stable page title and label")
            continue

        label = title.group(2)
        labels.add(label)
        if chr(96) in title.group(1):
            fail(errors, f"{source.relative_to(root)}: code markup is not allowed in page titles")
        if '<div class="distribution-meta"' not in text:
            fail(errors, f"{source.relative_to(root)}: missing distribution summary")
        if "## Definition and properties" not in text:
            fail(errors, f"{source.relative_to(root)}: missing definition section")
        if text.count(r"\f[") != text.count(r"\f]"):
            fail(errors, f"{source.relative_to(root)}: unbalanced display-math delimiters")
        if text.count(r"\f$") % 2:
            fail(errors, f"{source.relative_to(root)}: unbalanced inline-math delimiters")
        if re.search(r"(?<!\\)\bqquad\b", text):
            fail(errors, f"{source.relative_to(root)}: bare 'qquad' is missing a TeX backslash")

        rendered = html_dir / f"{label}.html"
        if not rendered.is_file():
            fail(errors, f"{rendered}: Doxygen did not generate the distribution page")
            continue
        output = rendered.read_text(encoding="utf-8")
        if "&lt;tt&gt;" in output or "&lt;/tt&gt;" in output:
            fail(errors, f"{rendered}: page title contains escaped code markup")
        if r"\f[" in output or r"\f]" in output or r"\f$" in output:
            fail(errors, f"{rendered}: raw Doxygen math delimiter leaked into HTML")
        displays_open = len(re.findall(r"(?m)^\s*\\\[\s*$", output))
        displays_close = len(re.findall(r"(?m)^\s*\\\]\s*$", output))
        if displays_open != displays_close:
            fail(errors, f"{rendered}: unbalanced MathJax display delimiters")
        if output.count(r"\(") != output.count(r"\)"):
            fail(errors, f"{rendered}: unbalanced MathJax inline delimiters")

    for index_name in ("continuous", "discrete", "bivariate", "other"):
        source = (root / f"docs/reference/{index_name}.md").read_text(encoding="utf-8")
        for target in re.findall(r'href="([a-z0-9_]+_distribution\.html)"', source):
            card_labels.add(target.removesuffix(".html"))
            if not (html_dir / target).is_file():
                fail(errors, f"docs/reference/{index_name}.md: broken card link to {target}")

    for label in sorted(labels - card_labels):
        fail(errors, f"{label}.html: distribution page is not linked from a family index")
    for label in sorted(card_labels - labels):
        fail(errors, f"{label}.html: family card does not have a distribution source page")

    for rendered in sorted(html_dir.glob("*_distribution.html")):
        output = rendered.read_text(encoding="utf-8")
        for asset in re.findall(r'<img[^>]+src="(images/[^"]+\.png)"', output):
            if not (html_dir / asset).is_file():
                fail(errors, f"{rendered}: missing plot asset {asset}")

    if errors:
        print("Rendered reference validation failed:", file=sys.stderr)
        for error in errors:
            print(f"  {error}", file=sys.stderr)
        return 1

    print(f"Validated structure, TeX delimiters, links, and assets for {len(labels)} distribution pages.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
