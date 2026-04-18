"""Top-level entrypoint for MCP hosts that load `server.py:mcp` as a file.

The canonical way to run the server is via the `oligomcp-mcp` console
script (see `oligomcp.mcp_server:main`), which uses the stdio transport
and is what Claude Code / Claude Desktop / Cursor / etc. expect. This
file exists as a fallback for any host that prefers a file:object
entrypoint — it just re-exports the FastMCP instance through the
installed package's absolute import path, so relative imports inside
`src/oligomcp/mcp_server.py` resolve cleanly.
"""
from oligomcp.mcp_server import mcp

__all__ = ["mcp"]
