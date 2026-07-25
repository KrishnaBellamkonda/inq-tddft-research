#!/usr/bin/env python3
"""Standalone tests for build_index.py (no pytest — keeps the skill shippable).

Run:  python3 test_build_index.py     # prints OK or raises AssertionError
"""
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import build_index as bi  # noqa: E402


def _campaign(area, slug, status, tasks_yaml, **extra):
    extra_lines = "".join(f"{k}: {v}\n" for k, v in extra.items())
    return (
        "---\n"
        f"id: {slug}\n"
        f"area: {area}\n"
        f"title: {slug.replace('-', ' ').title()}\n"
        f"status: {status}\n"
        f"{extra_lines}"
        "tasks:\n"
        f"{tasks_yaml}"
        "---\n\n# body\n"
    )


def build_fixture(root: Path):
    a = root / "alpha"
    b = root / "beta"
    a.mkdir()
    b.mkdir()

    # flow-style tasks, 2/3 done, running
    (a / "flow.md").write_text(_campaign(
        "alpha", "flow-camp", "running",
        '  - { name: "t1", done: true }\n'
        '  - { name: "t2", done: true }\n'
        '  - { name: "t3", done: false }\n',
        handover="docs/handovers/flow.md",
    ))
    # block-style tasks, 1/2 done, blocked w/ reason
    (b / "block.md").write_text(_campaign(
        "beta", "block-camp", "blocked",
        '  - name: "t1"\n    done: true\n'
        '  - name: "t2"\n    done: false\n',
        blocked_reason="waiting on GS",
    ))
    # done, 1/1
    (a / "done.md").write_text(_campaign(
        "alpha", "done-camp", "done",
        '  - { name: "t1", done: true }\n',
    ))
    # paused, 0/1, with a reason
    (b / "paused.md").write_text(_campaign(
        "beta", "paused-camp", "paused",
        '  - { name: "t1", done: false }\n',
        blocked_reason="parked till next week",
    ))
    # files that must be SKIPPED
    (root / "template.md").write_text("---\nfoo: bar\n---\n")     # template
    (root / "INDEX.md").write_text("stale\n")                     # the output
    (a / "notes.md").write_text("# just notes, no frontmatter\n")  # no id
    (a / "noid.md").write_text("---\narea: alpha\n---\n# no id\n")  # fm, no id


def main() -> int:
    with tempfile.TemporaryDirectory() as td:
        root = Path(td)
        build_fixture(root)
        camps = bi.collect(root)

        # Exactly the 4 real campaigns are picked up (4 non-campaigns skipped).
        ids = sorted(c["id"] for c in camps)
        assert ids == ["block-camp", "done-camp", "flow-camp", "paused-camp"], ids

        by_id = {c["id"]: c for c in camps}
        assert (by_id["flow-camp"]["_done"], by_id["flow-camp"]["_total"]) == (2, 3)
        assert (by_id["block-camp"]["_done"], by_id["block-camp"]["_total"]) == (1, 2)
        assert (by_id["done-camp"]["_done"], by_id["done-camp"]["_total"]) == (1, 1)

        out = bi.render(camps)

        # Portfolio header counts each status once.
        assert "1 running · 1 blocked · 1 paused · 1 done  (4 campaigns)" in out, out

        # Grouping order: running → blocked → paused → … → done.
        assert (out.index("🟢 Running") < out.index("⛔ Blocked")
                < out.index("⏸️ Paused") < out.index("✅ Done"))

        # Progress + blocked/paused reason rendered; handover link rewritten relative.
        assert "2/3" in out and "1/2 — *waiting on GS*" in out
        assert "0/1 — *parked till next week*" in out
        assert "[↗](../handovers/flow.md)" in out
        # done-camp has no handover -> em dash.
        assert "| done-camp" not in out  # title is shown, not id
        assert "Done Camp" in out

        # Skipped files never appear.
        assert "notes" not in out and "noid" not in out and "template" not in out

    print("OK — all build_index tests passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
