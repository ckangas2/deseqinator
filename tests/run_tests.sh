#!/usr/bin/env bash
# Unit tests: one planted failure per check. Seconds; no STAR, R or network.
set -euo pipefail
cd "$(dirname "$0")/.."
python3 -m pytest tests -q "$@"
