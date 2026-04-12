"""CLI entrypoint for OligoClaude."""
from __future__ import annotations

import argparse
import getpass
import sys
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="oligoclaude",
        description="Predict ASO efficacy via AlphaGenome and SpliceAI.",
    )
    sub = p.add_subparsers(dest="cmd", required=True)

    r = sub.add_parser("run", help="Run the full ASO prediction workflow.")
    r.add_argument("--config", type=Path, required=True)
    r.add_argument("--skip-alphagenome", action="store_true")
    r.add_argument("--skip-spliceai", action="store_true")
    r.add_argument(
        "--samples-max",
        type=int,
        default=20,
        help="Top/bottom ASOs to emit in compact BED (does not affect correlation).",
    )
    r.add_argument("--verbose", "-v", action="store_true")

    k = sub.add_parser("set-api-key", help="Save AlphaGenome API key (mode 0600).")
    k.add_argument("key", nargs="?", help="Key value; prompted if omitted.")

    sub.add_parser("clear-api-key", help="Remove the saved AlphaGenome API key.")

    for name, helptxt in [
        ("fetch-genome", "Download GRCh38 FASTA to ~/.oligoclaude/genomes/."),
        ("fetch-spliceai-weights", "Download the MANE-10000nt ensemble."),
    ]:
        f = sub.add_parser(name, help=helptxt)
        f.add_argument("--cache-dir", type=Path, default=None)

    return p


def _cmd_run(args) -> int:
    from .workflow import run_workflow

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


def _cmd_set_api_key(args) -> int:
    from .resources import save_alphagenome_api_key

    key = args.key or getpass.getpass("AlphaGenome API key (hidden): ")
    path = save_alphagenome_api_key(key)
    print(f"Saved AlphaGenome API key to {path} (mode 0600).")
    return 0


def _cmd_clear_api_key(args) -> int:
    from .resources import CRED_PATH, clear_alphagenome_api_key

    if clear_alphagenome_api_key():
        print(f"Removed AlphaGenome API key from {CRED_PATH}.")
    else:
        print("No saved AlphaGenome API key to remove.")
    return 0


def _cmd_fetch_genome(args) -> int:
    from .resources import ensure_hg38_fasta

    print(f"GRCh38 FASTA ready at: {ensure_hg38_fasta(args.cache_dir, verbose=True)}")
    return 0


def _cmd_fetch_spliceai(args) -> int:
    from .resources import ensure_spliceai_weights

    print(f"SpliceAI weights ready at: {ensure_spliceai_weights(args.cache_dir, verbose=True)}")
    return 0


_HANDLERS = {
    "run": _cmd_run,
    "set-api-key": _cmd_set_api_key,
    "clear-api-key": _cmd_clear_api_key,
    "fetch-genome": _cmd_fetch_genome,
    "fetch-spliceai-weights": _cmd_fetch_spliceai,
}


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    handler = _HANDLERS.get(args.cmd)
    if handler is None:
        parser.print_help()
        return 1
    return handler(args)


if __name__ == "__main__":
    sys.exit(main())
