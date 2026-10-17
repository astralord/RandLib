from __future__ import annotations

import unittest

from release_notes import extract


CHANGELOG = """\
# Changelog

## [2.0.0] - Unreleased

Release candidate notes.

## [2.0.0-beta.1] - 2026-10-17

Beta notes.

## [1.2.0] - 2026-10-16

Previous release notes.
"""


class ReleaseNotesTest(unittest.TestCase):
    def test_extracts_exact_release(self) -> None:
        self.assertEqual(extract(CHANGELOG, "v2.0.0"), "Release candidate notes.\n")

    def test_prefers_exact_prerelease(self) -> None:
        self.assertEqual(extract(CHANGELOG, "v2.0.0-beta.1"), "Beta notes.\n")

    def test_falls_back_to_base_release_for_release_candidate(self) -> None:
        self.assertEqual(extract(CHANGELOG, "v2.0.0-rc.1"), "Release candidate notes.\n")

    def test_rejects_unknown_release(self) -> None:
        with self.assertRaisesRegex(ValueError, r"no \[3\.0\.0-rc\.1\] section"):
            extract(CHANGELOG, "v3.0.0-rc.1")


if __name__ == "__main__":
    unittest.main()
