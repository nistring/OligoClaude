"""Snapshot the list of ontology terms AlphaGenome supports.

Queries `DnaClient.output_metadata(Organism.HOMO_SAPIENS)`, iterates over
every output type that ships biosample-level ontology metadata, and
emits a deduplicated TSV where each row is a unique `ontology_curie`
(e.g. ``CL:0000127``, ``UBERON:0002240``) annotated with the
human-readable biosample info and the set of AlphaGenome output types
that actually have at least one track labelled with that term.

The resulting TSV is committed to `data/alphagenome_ontology_terms.tsv`
so users can grep it offline when choosing `ontology_terms` values for
their config — the `predict_aso_efficacy*` tools accept these CURIEs
directly.

Usage:
    ALPHAGENOME_API_KEY=... python scripts/fetch_ontology_terms.py
    # or if the key is stored in ~/.oligomcp/credentials.json:
    python scripts/fetch_ontology_terms.py

Writes to:
    data/alphagenome_ontology_terms.tsv      (one row per CURIE)
    data/alphagenome_ontology_terms.meta.txt (generation timestamp)
"""
from __future__ import annotations

import datetime as _dt
import sys
from pathlib import Path

import pandas as pd


_OUTPUT_ATTRS = [
    "rna_seq",
    "cage",
    "procap",
    "atac",
    "dnase",
    "chip_histone",
    "chip_tf",
    "splice_sites",
    "splice_site_usage",
    "splice_junctions",
    "contact_maps",
]


def build_ontology_table() -> pd.DataFrame:
    """Return a deduplicated ontology-term table across all AlphaGenome outputs."""
    # Delay heavy imports so `-h` / `--help` doesn't pay for them.
    from alphagenome.models.dna_client import Organism, create

    from oligomcp.resources import require_alphagenome_api_key

    client = create(require_alphagenome_api_key())
    meta = client.output_metadata(Organism.HOMO_SAPIENS)

    frames: list[pd.DataFrame] = []
    for attr in _OUTPUT_ATTRS:
        df = getattr(meta, attr, None)
        if df is None or not hasattr(df, "columns"):
            continue
        if "ontology_curie" not in df.columns:
            continue
        keep = ["ontology_curie"]
        for optional in (
            "biosample_name",
            "biosample_type",
            "biosample_life_stage",
            "gtex_tissue",
        ):
            if optional in df.columns:
                keep.append(optional)
        piece = df[keep].dropna(subset=["ontology_curie"]).copy()
        piece["__output_type"] = attr
        piece["__track_count"] = 1  # to sum later
        frames.append(piece)

    if not frames:
        raise RuntimeError("No ontology metadata found — is the API key valid?")

    cat = pd.concat(frames, ignore_index=True)

    # Per CURIE: combine the biosample annotation (first non-null wins, since
    # the same CURIE gets identical labels across output types) and record
    # which output types actually have tracks for it.
    grouped = (
        cat.groupby("ontology_curie", as_index=False)
        .agg(
            biosample_name=("biosample_name", "first") if "biosample_name" in cat.columns else ("ontology_curie", "first"),
            biosample_type=("biosample_type", "first") if "biosample_type" in cat.columns else ("ontology_curie", "first"),
            biosample_life_stage=(
                "biosample_life_stage", "first"
            ) if "biosample_life_stage" in cat.columns else ("ontology_curie", "first"),
            gtex_tissue=("gtex_tissue", "first") if "gtex_tissue" in cat.columns else ("ontology_curie", "first"),
            available_outputs=(
                "__output_type",
                lambda s: ",".join(sorted(set(s))),
            ),
            total_tracks=("__track_count", "sum"),
        )
        .sort_values("ontology_curie")
        .reset_index(drop=True)
    )
    return grouped


def main() -> int:
    project_root = Path(__file__).resolve().parent.parent
    data_dir = project_root / "data"
    data_dir.mkdir(parents=True, exist_ok=True)

    tsv_path = data_dir / "alphagenome_ontology_terms.tsv"
    meta_path = data_dir / "alphagenome_ontology_terms.meta.txt"

    print("Fetching AlphaGenome ontology metadata …")
    df = build_ontology_table()
    df.to_csv(tsv_path, sep="\t", index=False)

    try:
        import alphagenome  # for version stamp

        ag_version = getattr(alphagenome, "__version__", "unknown")
    except Exception:
        ag_version = "unknown"

    meta_lines = [
        f"Generated at: {_dt.datetime.now(_dt.timezone.utc).isoformat()}",
        f"AlphaGenome version: {ag_version}",
        f"Unique ontology CURIEs: {len(df)}",
        f"Outputs covered: {sorted(set(','.join(df['available_outputs']).split(',')))}",
        "",
        "Usage: drop any CURIE (e.g. CL:0000127, UBERON:0002240) into",
        "the `ontology_terms` list of an OligoMCP config. Empty list =",
        "average over every track AlphaGenome returns.",
    ]
    meta_path.write_text("\n".join(meta_lines) + "\n", encoding="utf-8")

    print(f"Wrote {tsv_path} ({len(df)} unique CURIEs)")
    print(f"Wrote {meta_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
