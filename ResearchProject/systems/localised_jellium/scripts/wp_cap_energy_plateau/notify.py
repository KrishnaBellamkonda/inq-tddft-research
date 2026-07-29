#!/usr/bin/env python3
"""notify.py "subject" "body" [attachment ...] — thin email sender for the
autonomous orchestrator. Threads all messages under one Gmail family by reusing
a Message-ID stored in .mail_thread. Non-fatal: never raises to the caller.
"""
import sys
from pathlib import Path

FAMILY = "[wp-cap-energy-plateau]"
TO = "chiddukanna@gmail.com"
THREAD_FILE = Path(__file__).with_name(".mail_thread")


def main() -> int:
    if len(sys.argv) < 3:
        print("usage: notify.py subject body [attachment ...]")
        return 0
    subject = f"{FAMILY} {sys.argv[1]}"
    body = sys.argv[2]
    attachments = [a for a in sys.argv[3:] if Path(a).exists()]
    try:
        from inqview.email import send_run_email
        prior = THREAD_FILE.read_text().strip() if THREAD_FILE.exists() else None
        msg_id = send_run_email(
            subject, body, attachments=attachments or None, to=TO,
            in_reply_to=prior, references=[prior] if prior else None,
        )
        if not THREAD_FILE.exists() and msg_id:
            THREAD_FILE.write_text(msg_id)
        print(f"emailed: {subject}")
    except Exception as exc:  # never let email failure break the chain
        print(f"[notify] email failed ({exc}); continuing")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
