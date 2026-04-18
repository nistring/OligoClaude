"""Top-level entrypoint for Prefect Horizon / FastMCP Cloud deployment.

Horizon loads the configured entrypoint file directly as a script, so
relative imports inside `src/oligomcp/mcp_server.py` can't resolve.
This file just re-exports the `mcp` object via the installed package's
absolute import path, which works because Horizon installs dependencies
(including this package itself) from `pyproject.toml`.

Horizon configuration:
    Entrypoint: server.py:mcp
    Requirements: pyproject.toml (auto-detected)

Local stdio use is unchanged — the `oligomcp-mcp` console script still
calls `oligomcp.mcp_server:main`.
"""
from oligomcp.mcp_server import mcp

__all__ = ["mcp"]
