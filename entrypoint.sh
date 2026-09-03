#!/usr/bin/env bash
# Entrypoint for the API backend. Replaces bcai's "client.py & streamlit" combo
# with a single uvicorn process — the API IS the long-running service now.
set -e

# The bcai env is already on PATH (see Dockerfile), so `uvicorn` resolves to it.
cd /workspace

# Named Docker volumes mount with root ownership; create directories at runtime
mkdir -p /workspace/data
mkdir -p /workspace/user_data

# External MCP servers are bind-mounted from source (not baked into the image),
# so each needs a one-time npm build. Skipped once build/index.js exists.
# Loops over whatever's actually present under mcp_external/ instead of a
# hardcoded server list, so adding/removing a server needs no entrypoint change.
shopt -s nullglob
for dir in /workspace/mcp_external/*/; do
    server=$(basename "$dir")
    if [ -f "$dir/package.json" ] && [ ! -f "$dir/build/index.js" ]; then
        echo "Building $server..."
        (cd "$dir" && npm install && npm run build)
    fi
done

echo "Starting Drug Discovery Platform API on :8000"
exec uvicorn api.server:app --host 0.0.0.0 --port 8000
