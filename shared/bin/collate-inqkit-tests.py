#!/usr/bin/env python3
"""Collate inqkit Catch2 test results from ctest JUnit XML into a table.

Reads the per-tier `results.junit.xml` written by
`shared/bin/run-inqkit-tests.slurm` (ctest --output-junit) and prints a
per-test table plus per-tier and overall totals.

Parsing the XML rather than scraping ctest console text is deliberate: the
console summary reports only counts, and a test that never built at all is
invisible there. Here a tier whose XML is missing is reported explicitly as
NO RESULTS instead of silently contributing zero.

Usage:
    python3 shared/bin/collate-inqkit-tests.py [--md]

Pure stdlib, 3.6-compatible (bare `python3` on CSD3 is 3.6.8).
"""

import os
import sys
import xml.etree.ElementTree as ET

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))

TIERS = [
    ("pure", os.path.join(REPO, "inq-stack/tests/include/build")),
    ("engine", os.path.join(REPO, "inq-stack/tests/include/engine/build")),
]


def classify(tc):
    """Map a <testcase> to one of pass/fail/skip.

    ctest's JUnit writer marks a non-run test with status="notrun" and emits a
    <skipped/> child; a failure gets a <failure> child. Absence of both means it
    ran and passed.
    """
    if tc.find("failure") is not None or tc.find("error") is not None:
        return "FAIL"
    if tc.find("skipped") is not None or tc.get("status") == "notrun":
        return "SKIP"
    return "PASS"


def failure_reason(tc):
    node = tc.find("failure")
    if node is None:
        node = tc.find("error")
    if node is None:
        return ""
    msg = (node.get("message") or "").strip().replace("\n", " ")
    return msg[:90]


def load_tier(label, build_dir):
    xml = os.path.join(build_dir, "results.junit.xml")
    if not os.path.exists(xml):
        return None, "no results.junit.xml (build failed, or ctest never ran)"
    try:
        root = ET.parse(xml).getroot()
    except ET.ParseError as e:
        return None, "unparseable XML: {}".format(e)
    # ctest writes a single <testsuite>; tolerate a <testsuites> wrapper.
    suites = [root] if root.tag == "testsuite" else list(root.iter("testsuite"))
    rows = []
    for s in suites:
        for tc in s.iter("testcase"):
            try:
                t = float(tc.get("time") or 0.0)
            except ValueError:
                t = 0.0
            rows.append((tc.get("name") or "?", classify(tc), t, failure_reason(tc)))
    return rows, None


def main():
    md = "--md" in sys.argv
    all_rows = []
    problems = []

    for label, bdir in TIERS:
        rows, err = load_tier(label, bdir)
        if rows is None:
            problems.append((label, err))
            continue
        for name, status, t, reason in rows:
            all_rows.append((label, name, status, t, reason))

    if md:
        print("| Tier | Test | Result | Time (s) |")
        print("|---|---|---|---|")
        for label, name, status, t, _ in all_rows:
            print("| {} | `{}` | {} | {:.2f} |".format(label, name, status, t))
    else:
        print("{:<8} {:<44} {:<6} {:>8}".format("TIER", "TEST", "RESULT", "TIME(s)"))
        print("-" * 70)
        for label, name, status, t, reason in all_rows:
            print("{:<8} {:<44} {:<6} {:>8.2f}".format(label, name, status, t))
            if reason:
                print("{:<8} {:<44} {}".format("", "  ->", reason))

    print()
    grand = {"PASS": 0, "FAIL": 0, "SKIP": 0}
    for label, _ in TIERS:
        sub = [r for r in all_rows if r[0] == label]
        if not sub:
            continue
        c = {"PASS": 0, "FAIL": 0, "SKIP": 0}
        for r in sub:
            c[r[2]] += 1
            grand[r[2]] += 1
        print("{:<8} {:>3} tests: {:>3} pass, {:>3} fail, {:>3} skip   ({:.1f}s)".format(
            label, len(sub), c["PASS"], c["FAIL"], c["SKIP"], sum(r[3] for r in sub)))

    total = sum(grand.values())
    print()
    print("TOTAL    {:>3} tests: {:>3} pass, {:>3} fail, {:>3} skip".format(
        total, grand["PASS"], grand["FAIL"], grand["SKIP"]))

    for label, err in problems:
        print("MISSING  tier '{}': {}".format(label, err))

    failed = [r for r in all_rows if r[2] == "FAIL"]
    if failed:
        print()
        print("FAILED TESTS:")
        for label, name, _, _, reason in failed:
            print("  [{}] {}{}".format(label, name, "  -- " + reason if reason else ""))

    # Non-zero if anything failed OR a tier produced no results at all.
    return 1 if (failed or problems) else 0


if __name__ == "__main__":
    raise SystemExit(main())
