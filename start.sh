#!/usr/bin/env bash
# Container entrypoint: install dependencies at startup, then serve the API.
set -euo pipefail

cd "$(dirname "$0")"

echo "[start] installing dependencies from requirements.txt ..."
pip install --no-cache-dir -r requirements.txt

HOST="${HOST:-0.0.0.0}"
PORT="${PORT:-8000}"
echo "[start] launching uvicorn on ${HOST}:${PORT}"
exec uvicorn app.main:app --host "${HOST}" --port "${PORT}"
