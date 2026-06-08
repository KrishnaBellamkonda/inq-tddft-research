"""Send email via Gmail SMTP (App Password authentication).

Why this module exists
----------------------
The claude.ai Gmail MCP connector exposes ``create_draft`` only — no
``send_message``. For the 2026-05 jellium campaign we want analyse.py
post-processing pipelines to *send* emails programmatically, not just
queue drafts. This module uses Python's stdlib smtplib + a Google
App Password to send through Gmail's SMTPS endpoint.

One-time setup
--------------
1. Enable 2-Step Verification on your Google account (required to
   generate App Passwords): https://myaccount.google.com/security
2. Generate an App Password (16-character token, spaces optional):
   https://myaccount.google.com/apppasswords
3. Run::

       python -m inqview.email setup

   This prompts for your Gmail address and App Password and stores
   them at ``~/.config/inqview/gmail_credentials.json`` with
   ``chmod 600``. The credentials file lives outside the repo so it
   cannot be accidentally committed.

Usage
-----
::

    from inqview.email import send_run_email

    # Initial email of a thread.
    msg_id = send_run_email(
        subject="[jellium-knudsen-sweep] E=700 eV pair (1/5)",
        body="Pair summary follows.\\n\\nConfiguration table:\\n...",
        attachments=[
            "results/analysis/density_rt_delta_xz.gif",
            "results/analysis/energy_components.png",
            ...
        ],
    )

    # Reply (threaded under the same Gmail conversation).
    next_id = send_run_email(
        subject="[jellium-knudsen-sweep] E=800 eV pair (2/5)",
        body="...",
        attachments=[...],
        in_reply_to=msg_id,
        references=[msg_id],
    )

Returns
-------
Each call returns the Message-ID header of the sent email. Pass it as
``in_reply_to`` (and append to ``references``) to thread the next
email under the same Gmail conversation.

Smoke test
----------
::

    python -m inqview.email                       # sends test email to self

Limitations
-----------
- Gmail's SMTP limit is 500 recipients/day for free accounts. For the
  jellium campaign (one self-send per pair, ~30 pairs across the
  campaign) we are far below that.
- Attachment size cap is 25 MB per email (Gmail limit). At ~6 MB per
  300-frame GIF + a few PNGs, a pair email lands around 15-20 MB —
  tight but workable. For larger payloads, upload the GIF to a shared
  Drive folder and link to it in the body instead.
- The App Password is stored in plain text on disk at chmod 600.
  Anyone with shell access to your account can read it. This is the
  same threat model as ssh private keys; acceptable for a personal
  workstation, NOT for a shared/multi-user system.
"""

from __future__ import annotations

import email.utils
import json
import smtplib
import sys
from email import encoders
from email.mime.base import MIMEBase
from email.mime.image import MIMEImage
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from pathlib import Path


CREDENTIALS_PATH = Path.home() / ".config" / "inqview" / "gmail_credentials.json"

SMTP_HOST = "smtp.gmail.com"
SMTP_PORT_SSL = 465

_IMAGE_SUFFIXES = {".png", ".jpg", ".jpeg", ".gif"}


def _load_credentials() -> tuple[str, str]:
    if not CREDENTIALS_PATH.exists():
        raise FileNotFoundError(
            f"Gmail credentials not found at {CREDENTIALS_PATH}.\n"
            f"Run:  python -m inqview.email setup\n"
            f"(See module docstring for the one-time setup steps.)"
        )
    try:
        with open(CREDENTIALS_PATH) as f:
            creds = json.load(f)
    except json.JSONDecodeError as exc:
        raise ValueError(
            f"Credentials file at {CREDENTIALS_PATH} is not valid JSON: {exc}"
        ) from exc
    try:
        return creds["email"], creds["app_password"]
    except KeyError as exc:
        raise KeyError(
            f"Credentials file missing key {exc}; expected 'email' and "
            f"'app_password'. Re-run:  python -m inqview.email setup"
        ) from exc


def _attach_file(msg: MIMEMultipart, path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Attachment not found: {path}")

    suffix = path.suffix.lower()
    with open(path, "rb") as f:
        data = f.read()

    if suffix in _IMAGE_SUFFIXES:
        # MIMEImage auto-detects png/jpeg/gif from the bytes header.
        part = MIMEImage(data, name=path.name)
    else:
        part = MIMEBase("application", "octet-stream")
        part.set_payload(data)
        encoders.encode_base64(part)

    part.add_header(
        "Content-Disposition",
        f"attachment; filename={path.name}",
    )
    msg.attach(part)


def send_run_email(
    subject: str,
    body: str,
    attachments: list[str | Path] | None = None,
    to: str | None = None,
    in_reply_to: str | None = None,
    references: list[str] | None = None,
    html_body: str | None = None,
) -> str:
    """Send an email via Gmail SMTP and return its Message-ID.

    Parameters
    ----------
    subject
        Subject line. For threading within a Gmail conversation, keep the
        same family prefix (e.g. ``[jellium-knudsen-sweep]``) across all
        emails in the thread.
    body
        Plain-text body (markdown formatting is rendered as plain text
        in Gmail). For richer formatting pass ``html_body`` instead/also.
    attachments
        File paths to attach. PNG/JPEG/GIF use MIMEImage; everything
        else is octet-stream + base64.
    to
        Recipient email. Defaults to the sender (self-to-self).
    in_reply_to
        Message-ID of the email being replied to. Required for proper
        Gmail threading after the first email of a family.
    references
        Accumulated References header. Pass the list of all prior
        Message-IDs in the thread (or just ``[in_reply_to]`` for
        single-step threading).
    html_body
        Optional HTML body. If provided, the email is sent as
        multipart/alternative with the plain ``body`` as the
        text/plain part.

    Returns
    -------
    str
        Message-ID header of the sent email. Pass it to the next call
        as ``in_reply_to`` to thread under the same Gmail conversation.
    """
    sender, app_password = _load_credentials()
    recipient = to or sender

    if html_body:
        msg = MIMEMultipart("mixed")
        alternative = MIMEMultipart("alternative")
        alternative.attach(MIMEText(body, "plain"))
        alternative.attach(MIMEText(html_body, "html"))
        msg.attach(alternative)
    else:
        msg = MIMEMultipart("mixed")
        msg.attach(MIMEText(body, "plain"))

    msg["From"] = sender
    msg["To"] = recipient
    msg["Subject"] = subject
    msg["Date"] = email.utils.formatdate(localtime=True)
    msg_id = email.utils.make_msgid(domain="inqview.gmail")
    msg["Message-ID"] = msg_id

    if in_reply_to:
        msg["In-Reply-To"] = in_reply_to
        chain = list(references or [])
        if in_reply_to not in chain:
            chain.append(in_reply_to)
        msg["References"] = " ".join(chain)

    for path in attachments or []:
        _attach_file(msg, Path(path))

    with smtplib.SMTP_SSL(SMTP_HOST, SMTP_PORT_SSL) as server:
        server.login(sender, app_password)
        server.send_message(msg)

    return msg_id


def _setup_credentials() -> None:
    """Interactively prompt for credentials and save to disk (chmod 600)."""
    print("inqview.email — Gmail SMTP credentials setup")
    print()
    print("Prerequisites:")
    print("  1. 2-Step Verification enabled on your Google account.")
    print("     https://myaccount.google.com/security")
    print("  2. App Password generated (16-character token).")
    print("     https://myaccount.google.com/apppasswords")
    print()
    email_addr = input("Gmail address: ").strip()
    app_password = input("App password (16 chars, spaces are stripped): ").strip()
    app_password = app_password.replace(" ", "")
    if len(app_password) != 16:
        print(
            f"WARNING: App password is {len(app_password)} chars after "
            f"stripping spaces; Google App Passwords are normally 16 "
            f"chars. Saving anyway — verify by sending the smoke test."
        )

    CREDENTIALS_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(CREDENTIALS_PATH, "w") as f:
        json.dump(
            {"email": email_addr, "app_password": app_password}, f, indent=2,
        )
    CREDENTIALS_PATH.chmod(0o600)
    print()
    print(f"Saved to {CREDENTIALS_PATH} (chmod 600).")
    print()
    print("Verify with:  python -m inqview.email")


def _smoke_test() -> None:
    """Send a self-to-self test email and print the Message-ID.

    Catches SMTPAuthenticationError specifically and prints the ranked
    candidate root causes + remediation order — so future users hitting
    `535 5.7.8 Bad Credentials` get a runbook in the failure output
    instead of having to look it up.
    """
    print("inqview.email — smoke test")
    print(f"Reading credentials from {CREDENTIALS_PATH}...")
    sender, _ = _load_credentials()
    print(f"Sending test email to {sender}...")
    try:
        msg_id = send_run_email(
            subject="[inqview.email] smoke test",
            body=(
                "If you can read this in your inbox, inqview.email is "
                "working end-to-end via Gmail SMTPS + App Password.\n\n"
                f"Sent at: {email.utils.formatdate(localtime=True)}\n"
            ),
        )
    except smtplib.SMTPAuthenticationError as exc:
        print()
        print(f"FAIL: SMTP auth rejected (code={exc.smtp_code}).")
        print(
            f"Gmail says: "
            f"{exc.smtp_error.decode('utf-8', 'replace').strip()}"
        )
        if exc.smtp_code == 535:
            print()
            print("Ranked candidate root causes (most-likely first):")
            print("  1. 2-Step Verification is NOT actually enabled on the")
            print("     Google account. Check at:")
            print("       https://myaccount.google.com/security")
            print("     The App Password page sometimes shows even when")
            print("     2SV is off, but the password won't authenticate.")
            print("  2. Typo in App Password. Regenerate at:")
            print("       https://myaccount.google.com/apppasswords")
            print("  3. App Password generated under a different Google")
            print("     account (multi-account browser confusion).")
            print("     Verify the top-right avatar on the App Password")
            print("     page is the same address you set up here:")
            print(f"       {sender}")
            print("  4. You typed your regular Gmail password instead of")
            print("     the App Password (Google retired non-OAuth")
            print("     password access in 2022).")
            print("  5. Workspace/Enterprise admin needs to allow App")
            print("     Passwords for your account (only if this is a")
            print("     Workspace Gmail; rare for personal accounts).")
            print()
            print("Fix order: verify (1), check security alerts in your")
            print("Gmail inbox, regenerate (2/3/4) in one step, then:")
            print("   python -m inqview.email setup")
            print("   python -m inqview.email")
        raise SystemExit(1)
    except Exception as exc:
        print(f"FAIL: {type(exc).__name__}: {exc}")
        raise SystemExit(1)
    print(f"Sent. Message-ID: {msg_id}")
    print(
        "Check your inbox at the address above. If it landed in spam, "
        "mark as not spam so the next send is delivered cleanly."
    )


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "setup":
        _setup_credentials()
    else:
        _smoke_test()
