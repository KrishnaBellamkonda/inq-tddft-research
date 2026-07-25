---
name: email-notifications
description: Use whenever sending an email about a simulation run, campaign phase, parameter sweep, or any result the user asked to be emailed ("email me after each phase", "send the result", "notify me when done", per-phase Gmail in autonomous campaigns). Sends via inqview.email.send_run_email (Gmail SMTP / App Password — the claude.ai Gmail MCP only creates drafts, it cannot send). ENFORCES the mandatory four-part result-email structure (hypothesis reminder → what was done → what the plot shows → conclusion) with at least one plot attached. Invoke before composing ANY outbound email from this project.
---

# Email notifications (result / phase emails)

Every outbound email from this project that reports a result — a run, a campaign
phase, a sweep point — MUST follow the mandatory structure below and carry at
least one plot. This is a hard contract (user decision, 2026-06-27): a result
email without the four parts and a figure is not done.

## When to use

- The user asked to be emailed about a run/phase/sweep ("email me after each
  phase", "send the result", "notify me when it finishes").
- An **autonomous campaign** emits a per-phase notification (see the `campaigns`
  skill `<notebook_contract>` / autonomous mechanics).
- Any time you send an email from this project.

## How emails are sent (mechanism)

The claude.ai Gmail MCP exposes `create_draft` only — it **cannot send**. Send
programmatically via the project module:

```python
from inqview.email import send_run_email   # inq-stack/python/inqview/email.py

send_run_email(
    subject="[localised-jellium GS] H0 — base WP vs classical E_total(0)",
    body=BODY,                       # the four-part plain-text body (below)
    attachments=["…/H0_base_difference.png"],   # >= 1 plot, REQUIRED
    to="chiddukanna@gmail.com",
    html_body=HTML,                  # optional richer formatting
)
```

Run it with the project venv: `/local/data/public/skcb2/tddft/venv/bin/python3`.
Credentials come from the Gmail App Password (`_load_credentials()` in the module).
PNG/JPEG/GIF attach inline as `MIMEImage`; everything else as a file attachment.

## MANDATORY structure — every result email contains ALL FOUR + a plot

1. **Hypothesis** — a succinct, one-or-two-sentence reminder of the falsifiable
   hypothesis this phase tests (lift it verbatim/condensed from the campaign
   `<hypothesis_ladder>`). The reader must know *what question* the email answers.
2. **What was done** — the runs/method in 2–4 lines: system, key parameters,
   how many runs, what was measured. Enough to reconstruct intent, not a wall.
3. **What the plot shows** — name the attached figure and say what its axes /
   curves are and what to look at. The plot is **required** (≥ 1 `.png`).
4. **Conclusion** — the verdict from THIS phase's results: confirmed / refuted /
   partial, the key number(s) at 2–3 s.f., and (if a ladder) what it means for the
   next phase. Never "looks reasonable" — state the falsifiable outcome.

Keep it succinct and human-scannable. Numbers at 2–3 s.f. (number-rounding rule);
carry units. Label inferences "Inference:".

## Subject convention

`[<area> <campaign-short>] <phase id> — <one-line result>`
e.g. `[localised-jellium GS] H0 — base gap is artifact-dominated (−650 eV near)`.

## Body template

```
HYPOTHESIS
  <one–two line falsifiable statement this phase tests>

WHAT WAS DONE
  - <system + key params>
  - <N runs / sweep points; what was measured>

PLOT (attached: <filename>.png)
  <what the axes/curves are; what to look at>

CONCLUSION
  <confirmed / refuted / partial> — <key number(s), 2–3 s.f., units>.
  <one line: implication for the next phase, if part of a ladder.>
```

## Checklist before sending

- [ ] All four sections present (hypothesis, what-was-done, plot, conclusion).
- [ ] At least one plot attached (`.png`), and the body names it + says what it shows.
- [ ] Numbers at 2–3 s.f. with units; inferences labelled.
- [ ] Subject follows the convention.
- [ ] Sent with `inqview.email.send_run_email` via the project venv (NOT a draft).

## Notes

- For many figures, attach the 1–2 highlight plots and link a Drive folder for the
  rest (keeps the email light) — see `inqview.email` docstring.
- Forbidden-words rule (commit-messages) does NOT apply to emails — but keep emails
  professional and result-focused.
