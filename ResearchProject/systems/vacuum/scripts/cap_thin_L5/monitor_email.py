#!/usr/bin/env python3
"""Autonomous progress emailer for the cap_thin_L5 sweep.

Polls the run dirs and emails threaded progress updates to the user as runs
complete (subject family `[cap-thin-L5]`), then a final "compute complete" note.
Figure generation + the figures email are done by the main session on completion.

    PYTHONPATH=.../inq-stack/python python3 monitor_email.py &
"""
import sys, time, subprocess
from pathlib import Path
sys.path.insert(0, "/local/data/public/skcb2/tddft/inq-stack/python")
from inqview.email import send_run_email

TO = "chiddukanna@gmail.com"
SWEEP = Path("/local/data/public/skcb2/tddft/ResearchProject/systems/vacuum/cap_thin_L5")
TOTAL = 33
HA_TO_EV = 27.211386245988

def parse(p):
    out = {}
    for ln in p.read_text().splitlines():
        k, _, v = ln.partition(' ')
        try: out[k] = float(v)
        except ValueError: out[k] = v
    return out

def done_recs():
    recs = []
    for d in sorted(SWEEP.glob('run_cap_*')):
        f = d / 'results/epsilon.txt'
        if f.exists():
            try: recs.append(parse(f))
            except Exception: pass
    return recs

def snapshot(recs):
    """Best ε so far per η, as text."""
    by_eta = {}
    for r in recs:
        e = round(r.get('eta_Ha', 0), 3)
        by_eta.setdefault(e, []).append(r)
    lines = []
    for e in sorted(by_eta):
        c = by_eta[e]
        best = min(c, key=lambda r: r.get('epsilon', 9))
        lines.append(f"  eta={e:+.2f} Ha: {len(c)} pts, best eps={best.get('epsilon',0):.3e} "
                     f"at E={best.get('E_eV',0):.1f} eV")
    return "\n".join(lines) if lines else "  (no runs parsed yet)"

def dispatch_alive():
    return subprocess.run(["pgrep", "-f", "cap_thin_L5/dispatch.py"],
                          capture_output=True).returncode == 0

def main():
    msgid = None
    refs = []
    last_sent = -1
    # milestones to email at (count of completed runs)
    milestones = [11, 18, 25, 32]
    sent = set()
    while True:
        recs = done_recs()
        n = len(recs)
        alive = dispatch_alive()
        # pick the highest milestone reached but not yet sent
        due = [m for m in milestones if n >= m and m not in sent]
        if due and n != last_sent:
            m = max(due)
            for x in [x for x in milestones if x <= m]:
                sent.add(x)
            subj = f"[cap-thin-L5] sweep progress {n}/{TOTAL} runs done"
            body = (f"Thin in-built CAP (L=5 Bohr) reflectivity sweep — progress update.\n\n"
                    f"Completed: {n}/{TOTAL} runs.\n\n"
                    f"Best reflection error so far, per depth:\n{snapshot(recs)}\n\n"
                    f"Lower eps = better absorption. Full ε(E) curves + density GIF follow "
                    f"when the sweep finishes.\nResults remain provisional pending the "
                    f"inq-study engine regression (Task #7).")
            try:
                mid = send_run_email(subj, body, to=TO, in_reply_to=msgid,
                                     references=refs or None)
                if msgid is None: msgid = mid
                refs.append(mid)
            except Exception as ex:
                print("email error:", ex, flush=True)
            last_sent = n
        if not alive and n >= TOTAL:
            subj = f"[cap-thin-L5] compute complete — {n}/{TOTAL} runs, building figures"
            body = (f"All {n}/{TOTAL} thin-CAP runs finished.\n\n"
                    f"Best reflection error per depth:\n{snapshot(recs)}\n\n"
                    f"Now assembling the ε(E) reflectivity curves and density GIF — "
                    f"those arrive in a follow-up email shortly.")
            try:
                send_run_email(subj, body, to=TO, in_reply_to=msgid, references=refs or None)
            except Exception as ex:
                print("final email error:", ex, flush=True)
            print(f"monitor done: {n}/{TOTAL}", flush=True)
            return
        if not alive and n < TOTAL:
            # dispatch died early; report and exit
            try:
                send_run_email(f"[cap-thin-L5] sweep STOPPED at {n}/{TOTAL}",
                               f"The dispatcher is no longer running but only {n}/{TOTAL} runs "
                               f"completed. Manual check needed.\n\n{snapshot(recs)}",
                               to=TO, in_reply_to=msgid, references=refs or None)
            except Exception as ex:
                print("stop email error:", ex, flush=True)
            print(f"monitor: dispatch died at {n}/{TOTAL}", flush=True)
            return
        time.sleep(45)

if __name__ == "__main__":
    main()
