#!/usr/bin/env bash
# Scaffold a dated submission folder for Chase's Journal.
# Usage: bash scripts/new_submission.sh "<short-slug>"
# Prints the path of the created submission folder.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SLUG="${1:-untitled}"
SLUG="$(echo "$SLUG" | tr '[:upper:] ' '[:lower:]-' | tr -cd 'a-z0-9-')"

MONTH="$(date +%Y-%m)"
DATE="$(date +%Y-%m-%d)"
MONTH_DIR="$ROOT/submissions/$MONTH"
mkdir -p "$MONTH_DIR"

# Index within the month = count of existing submission folders this month + 1
N=$(( $(find "$MONTH_DIR" -maxdepth 1 -mindepth 1 -type d 2>/dev/null | wc -l | tr -d ' ') + 1 ))

DIR="$MONTH_DIR/m${N}-${SLUG}"
mkdir -p "$DIR/figs"
cp "$ROOT/templates/submission.md" "$DIR/note.md"

cat > "$DIR/meta.json" <<JSON
{
  "title": "",
  "author_persona": "",
  "topic": "",
  "kind": "theoretical | computational | mixed",
  "contribution_type": "",
  "date": "$DATE",
  "status": "submitted",
  "resubmission_of": null
}
JSON

echo "$DIR"
