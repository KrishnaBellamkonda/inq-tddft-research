---
name: journal-writing
description: Use when appending a new entry to one of the project journals (`docs/journals/<journal>/`). Each journal lives in its own folder with `index.md` (table of contents) and one markdown file per entry. Captures one run with its config snapshot from `run_summary.txt`, the user's observation text, and any attached figures. Auto-copies images into `docs/journals/<journal>/attachments/<slug>/`.
---

# Journal Writing

This skill appends one entry per invocation to a project journal under
`docs/journals/<journal>/`. The journal is the user's permanent record of
what was run, when, with what config, and what they observed. The user
supplies the **observation text** — never invent it.

## Layout (current)

```
docs/journals/
├── tutorial/
│   ├── index.md                       (main page: TOC table of all entries)
│   ├── 2026-05-04_<slug>.md           (one file per entry)
│   ├── 2026-05-05_<slug>.md
│   └── attachments/
│       ├── 2026-05-04_<slug>/
│       │   └── plot1.png
│       └── 2026-05-05_<slug>/
│           └── ...
├── quantumkickextension/
│   ├── index.md
│   ├── 2026-05-03_<slug>.md
│   └── attachments/<slug>/...
└── researchproject/
    ├── index.md
    ├── 2026-05-03_<slug>.md
    └── attachments/<slug>/...
```

Each journal folder is self-contained: its own `index.md`, its own per-
entry markdown files, and its own `attachments/` subfolder. **All entries
within one journal share that journal's `attachments/` folder.** Entries
across journals do not share attachments.

Topical entries (notes, conclusions, cross-cutting writeups that are not
tied to a single run) live in the same folder, with a slug like
`<topic_name>` instead of `<date>_<run_name>`. They follow the same
template — the `Run path` / `Config snapshot` sections may be omitted or
replaced with a `Linked entries` list.

## When to use

Invoke this skill when:

- A run has finished (or has been declared complete enough to record) and
  the user has supplied an observation paragraph.
- The user asks to "write a journal entry", "add to the journal", or
  "record this run".
- The user asks to write a topical / cross-cutting entry that links to
  several earlier entries (treat as an entry without a `run_path`).

Do **not** invoke this skill to:

- Edit an existing entry silently (use `--update <slug>` only when the
  user explicitly says so — currently deferred / not implemented).
- Generate observation text on the user's behalf. Ask the user for it.

## Inputs requested from the user

When invoked, confirm or ask the user for:

1. **Journal name** — one of `quantumkickextension`, `researchproject`,
   `tutorial`. Pick by directory of the run path if not stated:
   - `QuantumKickExtension/...` → `quantumkickextension`
   - `ResearchProject/...` → `researchproject`
   - `Tutorial/...` → `tutorial`
2. **Run path** — absolute path to the run directory (the parent of
   `results/run_summary.txt`). For topical entries with no single run,
   skip this and accept a `Linked entries` list instead.
3. **Title** — short human-readable title for the entry (used in the
   index.md table; defaults to the `run` field of `run_summary.txt`).
4. **Observations text** — markdown blob in the user's voice. May
   include `![alt](path/to/figure.png)` references; absolute paths or
   paths relative to `<run_path>/results/` are both fine.
5. **Open questions / next steps** (optional) — markdown bullet list.
6. **Status** — one of `running`, `complete`, `failed` (defaults to
   `complete`).
7. **Slug override** (optional) — defaults to `YYYY-MM-DD_<run_name>`
   for run entries, or `<topic_name>` for topical entries.

## What this skill does (procedure)

### Step 1 — Collect inputs and verify paths

- Verify `<run_path>/results/run_summary.txt` exists (skip for topical
  entries). If a run was declared but the summary is missing, abort and
  ask the user.
- Determine `<run_name>` from `basename(<run_path>)` (or from the
  `run = ...` line in `run_summary.txt` if present).
- Determine `<slug>` from the user override or from
  `YYYY-MM-DD_<run_name>` (today's date, safe-shell version of name).
- Verify `docs/journals/<journal>/` exists; if not, create it with an
  empty `attachments/` and a fresh `index.md` (see Step 5 for template).

### Step 2 — Parse run_summary.txt into a markdown table

(Skip for topical entries.) Each line of `run_summary.txt` is
`key = value`. Build a two-column markdown table verbatim — no editorial
reordering:

```md
| Field | Value |
|---|---|
| run | run_propagate_v0p450_extensive |
| system | li_54_atom_bcc_supercell |
| ... | ... |
```

If a value contains a `|` or newline, escape or wrap in backticks.

### Step 3 — Resolve attachments

For every `![alt](src)` in the observation text:

- If `src` is already `attachments/<slug>/...` (relative under the
  journal folder), leave it alone.
- Otherwise, resolve `src`:
  - If absolute, use as-is.
  - If relative, resolve against `<run_path>/results/`.
  Then copy the file to
  `docs/journals/<journal>/attachments/<slug>/<basename(src)>` and
  rewrite the markdown reference (in the entry file) to
  `attachments/<slug>/<basename(src)>` (relative to the entry's own
  directory, so links work inside the journal folder).
- If the source file does not exist, abort with a clear error.

Use `cp` (not move) so the live `results/` tree is untouched.

### Step 4 — Write the per-entry markdown file

Write `docs/journals/<journal>/<slug>.md` (refusing to overwrite if it
already exists) using this template:

```md
# <slug>

**Title:** <user-supplied title>
**Run path:** `<absolute run_path>`              (omit for topical entries)
**Linked results:** `<absolute results path>`    (omit for topical entries)
**Status:** <status>

## Config snapshot

<two-column table from Step 2 — omit for topical entries>

## Observations

<the user's markdown blob, with image refs rewritten to attachments/<slug>/...>

## Open questions / next steps

<the user's bullet list — omit this section entirely if none supplied>
```

### Step 5 — Update the journal's `index.md`

Read or initialise `docs/journals/<journal>/index.md`. The file is a
single page with a level-1 heading and one TOC table (newest at the top).

If `index.md` does not exist, create it with this header:

```md
# Journal: <Journal display name>

Rolling table of contents for the `<journal>` journal. Each entry below
links to its own markdown file in this folder. Attachments live under
`attachments/<slug>/`.

| Date | Title | Slug | Status |
|---|---|---|---|
```

Append (insert at the top of the table, just after the header row) a
new row:

```md
| YYYY-MM-DD | <title> | [`<slug>`](<slug>.md) | <status> |
```

### Step 6 — Confirm

Print one line summary: the journal folder, the slug, the relative paths
of the new entry file and any attachments copied.

## Constraints

- **Append only.** Never edit an existing entry file. If
  `<slug>.md` already exists, refuse and surface the conflict to the
  user. Updates are deferred until an explicit `--update <slug>` flag is
  added later.
- **Never invent observation text or interpretation.** If the user has
  not supplied observations, ask. Do not summarise the data on their
  behalf.
- **Never write into the live `results/` tree.** All copies go to
  `docs/journals/<journal>/attachments/<slug>/`.
- **Absolute paths for run-related references.** `Run path` and
  `Linked results` use absolute paths.
- **One TOC table.** `index.md` has exactly one TOC table — no per-year
  sections, no narrative padding.
