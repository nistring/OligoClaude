"""Local MCP server exposing OligoClaude as a Claude connector.

Register via:
    pip install -e .
    claude mcp add oligoclaude -- oligoclaude-mcp

Or add to Claude Desktop ~/.config/Claude/claude_desktop_config.json:
    {
      "mcpServers": {
        "oligoclaude": {
          "command": "oligoclaude-mcp",
          "env": { "ALPHAGENOME_API_KEY": "your-key-here" }
        }
      }
    }

Then ask in any Claude session:
    "Predict ASO efficacy for /path/to/SETD5_e1.json"

Credential resolution order on startup:
    1. $ALPHAGENOME_API_KEY environment variable
    2. ~/.oligoclaude/credentials.json (saved via `oligoclaude set-api-key`)
    3. Legacy `dna_api_key` field inside the individual config JSON (discouraged)

The server warns on startup if no key is configured; callers that pass
`skip_alphagenome=true` can still use the SpliceAI-only path without a key.
"""
from __future__ import annotations

import sys
from pathlib import Path
from typing import Optional

from .resources import ENV_VAR, get_alphagenome_api_key
from .workflow import WorkflowResult, run_workflow


def _build_mcp():
    try:
        from mcp.server.fastmcp import FastMCP
    except Exception as e:
        raise RuntimeError(
            "The `mcp` package is not installed. Run `pip install mcp`."
        ) from e

    mcp = FastMCP("oligoclaude")

    @mcp.tool()
    def predict_aso_efficacy(
        config_path: str,
        skip_alphagenome: bool = False,
        skip_spliceai: bool = False,
        samples_max: int = 20,
    ) -> dict:
        """Predict ASO efficacy using AlphaGenome and SpliceAI.

        Runs the full OligoClaude workflow from a config JSON file:
        enumerates ASO candidates (sliding window or from experimental data),
        scores them with AlphaGenome and SpliceAI, writes per-source BED tracks
        for UCSC upload, and generates a correlation plot if experimental data
        is provided.

        Args:
            config_path: Absolute path to the OligoClaude JSON config file.
            skip_alphagenome: If True, skip AlphaGenome scoring (no API key needed).
            skip_spliceai: If True, skip SpliceAI scoring.
            samples_max: Number of top/bottom ASOs per BED track (default 20).

        Returns:
            A dict with these keys:
              scores_csv:        Path to the per-candidate scores CSV.
              bed_files:         List of generated BED file paths (UCSC custom tracks).
              correlation_plot:  Path to the PNG plot (if experimental data present).
              stats:             Pearson/Spearman correlations per source per exon.
              ucsc_instructions: Text instructions for uploading BEDs to UCSC.
              n_candidates:      Number of ASO candidates scored.
        """
        cfg_path = Path(config_path).expanduser().resolve()
        if not cfg_path.exists():
            raise FileNotFoundError(f"Config file not found: {cfg_path}")

        result: WorkflowResult = run_workflow(
            cfg_path,
            skip_alphagenome=skip_alphagenome,
            skip_spliceai=skip_spliceai,
            samples_max=samples_max,
        )
        return {
            "scores_csv": str(result.scores_csv) if result.scores_csv else None,
            "bed_files": [str(p) for p in result.bed_files],
            "correlation_plot": (
                str(result.correlation_plot) if result.correlation_plot else None
            ),
            "stats": _stats_to_json(result.stats),
            "ucsc_instructions": result.ucsc_instructions,
            "n_candidates": result.n_candidates,
        }

    return mcp


def _stats_to_json(stats: Optional[dict]) -> Optional[dict]:
    """Convert nested stats dict (tuples) to JSON-serializable nested dicts."""
    if not stats:
        return None
    out: dict[str, dict] = {}
    for exon, by_src in stats.items():
        out[exon] = {}
        for src, vals in by_src.items():
            r, p, rho, rp = vals
            out[exon][src] = {
                "pearson_r": r,
                "pearson_p": p,
                "spearman_rho": rho,
                "spearman_p": rp,
            }
    return out


def _check_startup_credentials() -> None:
    """Print a one-line advisory to stderr if no AlphaGenome key is configured.

    Does not abort — callers may still invoke `predict_aso_efficacy` with
    `skip_alphagenome=true` for SpliceAI-only runs.
    """
    if get_alphagenome_api_key() is None:
        sys.stderr.write(
            "[oligoclaude-mcp] WARNING: no AlphaGenome API key found.\n"
            f"  Set ${ENV_VAR} in the MCP server `env` block, or run "
            "`oligoclaude set-api-key <KEY>` once.\n"
            "  SpliceAI-only calls (skip_alphagenome=true) will still work.\n"
        )


def main() -> None:
    _check_startup_credentials()
    mcp = _build_mcp()
    mcp.run(transport="stdio")


if __name__ == "__main__":
    main()
