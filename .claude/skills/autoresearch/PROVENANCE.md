# Provenance

Installed 2026-07-13 from https://github.com/drivelineresearch/autoresearch-claude-code
(v1.1.0, MIT License — see LICENSE). Maintainer: Driveline Research. 327 stars at install.

Only the skill itself (SKILL.md + scripts/ar-log.sh) is installed — the repo's plugin
hooks (Stop/PreCompact/SessionStart) are NOT installed, so loop continuation is driven
by the session agent, not mechanically enforced. The JSONL/worklog/dashboard protocol is
followed as written so a future session with the full plugin can resume seamlessly.

First use: localised-jellium CAP energy-artifact fix campaign
(docs/campaigns/localised_jellium/cap-fix-experimentation.md).
