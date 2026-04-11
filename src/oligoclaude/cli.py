"""CLI entrypoint for OligoClaude."""
from __future__ import annotations

import argparse
import getpass
import sys
from pathlib import Path

from .workflow import run_workflow


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="oligoclaude",
        description="Predict ASO efficacy via AlphaGenome and SpliceAI, and "
        "compare to experimental data.",
    )
    sub = parser.add_subparsers(dest="cmd", required=True)

    p_run = sub.add_parser(
        "run", help="Run the full ASO prediction workflow for one config."
    )
    p_run.add_argument("--config", type=Path, required=True, help="Path to JSON config.")
    p_run.add_argument(
        "--skip-alphagenome", action="store_true", help="Skip AlphaGenome scoring."
    )
    p_run.add_argument(
        "--skip-spliceai", action="store_true", help="Skip SpliceAI scoring."
    )
    p_run.add_argument(
        "--samples-max",
        type=int,
        default=20,
        help="Top/bottom ASOs to emit in compact BED (does not affect correlation).",
    )
    p_run.add_argument("--verbose", "-v", action="store_true")

    p_key = sub.add_parser(
        "set-api-key",
        help="Save the AlphaGenome API key to ~/.oligoclaude/credentials.json (mode 0600).",
    )
    p_key.add_argument(
        "key",
        nargs="?",
        help="API key value. If omitted, you will be prompted (hidden input).",
    )

    sub.add_parser(
        "clear-api-key",
        help="Remove the saved AlphaGenome API key from ~/.oligoclaude/credentials.json.",
    )

    p_fg = sub.add_parser(
        "fetch-genome",
        help="Download GRCh38 primary assembly FASTA to ~/.oligoclaude/genomes/ (one-time).",
    )
    p_fg.add_argument(
        "--cache-dir",
        type=Path,
        default=None,
        help="Override the default cache directory.",
    )

    p_fs = sub.add_parser(
        "fetch-spliceai-weights",
        help="Download the OpenSpliceAI MANE-10000nt 5-model ensemble weights.",
    )
    p_fs.add_argument(
        "--cache-dir",
        type=Path,
        default=None,
        help="Override the default cache directory.",
    )

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.cmd == "run":
        result = run_workflow(
            config_path=args.config,
            skip_alphagenome=args.skip_alphagenome,
            skip_spliceai=args.skip_spliceai,
            samples_max=args.samples_max,
            verbose=args.verbose,
        )
        print()
        print(f"Candidates scored: {result.n_candidates}")
        print(f"Scores CSV: {result.scores_csv}")
        if result.correlation_plot:
            print(f"Correlation plot: {result.correlation_plot}")
        print()
        print(result.ucsc_instructions)
        return 0

    if args.cmd == "set-api-key":
        from .credentials import save_alphagenome_api_key

        key = args.key or getpass.getpass("AlphaGenome API key (hidden): ")
        path = save_alphagenome_api_key(key)
        print(f"Saved AlphaGenome API key to {path} (mode 0600).")
        return 0

    if args.cmd == "clear-api-key":
        from .credentials import CRED_PATH, clear_alphagenome_api_key

        if clear_alphagenome_api_key():
            print(f"Removed AlphaGenome API key from {CRED_PATH}.")
        else:
            print("No saved AlphaGenome API key to remove.")
        return 0

    if args.cmd == "fetch-genome":
        from .genome_fetch import ensure_hg38_fasta

        path = ensure_hg38_fasta(args.cache_dir, verbose=True)
        print(f"GRCh38 FASTA ready at: {path}")
        return 0

    if args.cmd == "fetch-spliceai-weights":
        from .spliceai_fetch import ensure_spliceai_weights

        path = ensure_spliceai_weights(args.cache_dir, verbose=True)
        print(f"SpliceAI weights ready at: {path}")
        return 0

    parser.print_help()
    return 1


if __name__ == "__main__":
    sys.exit(main())
