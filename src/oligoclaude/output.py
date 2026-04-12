"""Output writers: BED custom tracks for UCSC + per-exon correlation plots.

Both writers consume the same `candidates` + `scores` payload produced by
`predict.py`. BED export produces a compact (top/bottom 20) file and a full
file per score source; the correlation plot overlays one regression line
per score source per exon panel, plus an "all-exons" panel when multiple
exons are present. `open_ucsc_browser` auto-uploads compact BED tracks to
the UCSC Genome Browser via `hgct_customText` URL parameter.
"""
from __future__ import annotations

import urllib.parse
import webbrowser
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from .config import OligoConfig
from .core import AsoCandidate

_BED_COLS = [
    "chrom", "chromStart", "chromEnd", "name", "score",
    "strand", "thickStart", "thickEnd", "itemRgb",
]


# ---------- BED export ----------

def _build_rows(
    selected: list[tuple[AsoCandidate, float]],
    chrom: str,
    bed_strand: str,
    anchor: int,
    invert: bool,
) -> pd.DataFrame:
    """Render (candidate, score) pairs into a 9-column BED DataFrame."""
    if not selected:
        return pd.DataFrame(columns=_BED_COLS)

    cands, scores = zip(*selected)
    scores_arr = np.asarray(scores, dtype=np.float32)
    vals = np.clip((1 + scores_arr if invert else 1 - scores_arr) * 255, 1, 255).astype(int)
    if invert:
        item_rgb = [f"{int(v)},{int(v)},255" for v in vals]
        label = "bottom"
    else:
        item_rgb = [f"255,{int(v)},{int(v)}" for v in vals]
        label = "top"

    starts = [anchor + c.position for c in cands]
    ends = [anchor + c.position + c.length for c in cands]
    names = [f"{label}{i + 1}({s * 100:.2f}%)" for i, s in enumerate(scores_arr)]

    return pd.DataFrame(
        {
            "chrom": [chrom] * len(cands),
            "chromStart": starts,
            "chromEnd": ends,
            "name": names,
            "score": scores_arr,
            "strand": [bed_strand] * len(cands),
            "thickStart": starts,
            "thickEnd": ends,
            "itemRgb": item_rgb,
        }
    )


def _safe_concat(parts: list[pd.DataFrame]) -> pd.DataFrame:
    non_empty = [p for p in parts if len(p)]
    return pd.concat(non_empty, ignore_index=True) if non_empty else pd.DataFrame(columns=_BED_COLS)


def _write_track(path: Path, header: str, df: pd.DataFrame) -> None:
    with open(path, "w", newline="") as f:
        f.write(header + "\n")
        df.to_csv(f, sep="\t", index=False, header=False, lineterminator="\n")


def write_bed(
    results_dir: Path,
    config_name: str,
    source: str,
    chrom: str,
    strand: str,
    candidates: list[AsoCandidate],
    scores: np.ndarray,
    variant_interval_start: int,
    samples_max: int = 20,
) -> tuple[Path, Path]:
    """Write compact (top/bottom-N) and full BED files for one score source.

    `candidates[i].position` is the ref_seq offset; absolute genomic starts
    are `variant_interval_start + position`. The BED strand is inverted
    relative to the config strand (ASO binds the opposite strand). Positive
    scores get a red-gradient itemRgb, negative scores get blue.
    """
    results_dir = Path(results_dir)
    results_dir.mkdir(parents=True, exist_ok=True)
    bed_strand = "+" if strand == "-" else "-"

    pairs = list(zip(candidates, scores.tolist()))
    pos = sorted([p for p in pairs if p[1] > 0], key=lambda x: -x[1])
    neg = sorted([p for p in pairs if p[1] < 0], key=lambda x: x[1])

    def _render(pos_sel, neg_sel):
        return _safe_concat([
            _build_rows(pos_sel, chrom, bed_strand, variant_interval_start, invert=False),
            _build_rows(neg_sel, chrom, bed_strand, variant_interval_start, invert=True),
        ])

    compact_df = _render(pos[:samples_max], neg[:samples_max])
    full_df = _render(pos, neg)

    compact_path = results_dir / f"{config_name}_ASO_{source}.bed"
    full_path = results_dir / f"{config_name}_ASO_{source}_full.bed"
    _write_track(
        compact_path,
        f'track name={config_name}_ASO_{source} '
        f'description="ASO scores for {source}" visibility="pack" useScore=1',
        compact_df,
    )
    _write_track(
        full_path,
        f'track name={config_name}_ASO_{source}_full '
        f'description="Full ASO scores for {source}" visibility="pack" useScore=1',
        full_df,
    )
    return compact_path, full_path


def export_all(
    cfg: OligoConfig,
    chrom: str,
    variant_interval_start: int,
    candidates: list[AsoCandidate],
    all_scores: dict[str, np.ndarray],
    samples_max: int = 20,
) -> list[Path]:
    """Write a pair of BED files per score source."""
    written: list[Path] = []
    for source, arr in all_scores.items():
        if arr is None or len(arr) == 0:
            continue
        compact, full = write_bed(
            results_dir=cfg.results_dir,
            config_name=cfg.config_name,
            source=source,
            chrom=chrom,
            strand=cfg.strand,
            candidates=candidates,
            scores=np.asarray(arr, dtype=np.float32),
            variant_interval_start=variant_interval_start,
            samples_max=samples_max,
        )
        written.extend([compact, full])
    return written


def write_experimental_bed(
    results_dir: Path,
    config_name: str,
    chrom: str,
    strand: str,
    candidates: list[AsoCandidate],
    variant_interval_start: int,
) -> Optional[Path]:
    """Write a BED track from experimental measured (RT-PCR) values.

    Creates a green-gradient track where darker green = higher measured
    value. Returns None if no candidates have measured data.
    """
    exp = [(c, c.measured) for c in candidates if c.measured is not None]
    if not exp:
        return None

    results_dir = Path(results_dir)
    results_dir.mkdir(parents=True, exist_ok=True)
    bed_strand = "+" if strand == "-" else "-"

    vals = np.array([m for _, m in exp], dtype=np.float32)
    lo, hi = vals.min(), vals.max()
    norm = (vals - lo) / (hi - lo + 1e-12)
    intensities = np.clip((1 - norm) * 200 + 55, 55, 255).astype(int)

    rows = []
    for i, ((c, m), intensity) in enumerate(zip(exp, intensities)):
        s = variant_interval_start + c.position
        e = s + c.length
        rows.append({
            "chrom": chrom,
            "chromStart": s,
            "chromEnd": e,
            "name": f"{c.aso_id}({m:.1f})",
            "score": int(np.clip(m * 100, 0, 1000)),
            "strand": bed_strand,
            "thickStart": s,
            "thickEnd": e,
            "itemRgb": f"0,{int(intensity)},0",
        })

    path = results_dir / f"{config_name}_ASO_Measured.bed"
    _write_track(
        path,
        f'track name={config_name}_ASO_Measured '
        f'description="Measured RT-PCR" visibility="pack" useScore=1 '
        f'itemRgb="On"',
        pd.DataFrame(rows),
    )
    return path


def open_ucsc_browser(
    assembly: str,
    chrom: str,
    exon_start: int,
    exon_end: int,
    flank: tuple[int, int],
    bed_files: list[Path],
) -> Optional[str]:
    """Auto-open the UCSC Genome Browser with compact BED tracks loaded.

    Reads the compact (non-full) BED files, URL-encodes their content into
    the `hgct_customText` parameter, and opens the browser at the target
    exon ± flank. Returns the URL, or None if there are no tracks.
    """
    compact_files = [f for f in bed_files if not f.name.endswith("_full.bed")]
    if not compact_files:
        return None

    bed_content = "\n".join(f.read_text().rstrip() for f in compact_files if f.exists())
    if not bed_content.strip():
        return None

    position = f"{chrom}:{exon_start - flank[0]}-{exon_end + flank[1]}"
    params = urllib.parse.urlencode({
        "db": assembly,
        "position": position,
        "hgct_customText": bed_content,
    })
    url = f"https://genome.ucsc.edu/cgi-bin/hgTracks?{params}"

    try:
        webbrowser.open(url)
    except Exception:
        pass
    return url


# ---------- Correlation plot ----------

def _safe_corr(x: pd.Series, y: pd.Series):
    """Return (pearson_r, pearson_p, spearman_rho, spearman_p) with NaN on failure."""
    from scipy import stats

    mask = x.notna() & y.notna()
    if mask.sum() < 3:
        return (np.nan, np.nan, np.nan, np.nan)
    xc, yc = x[mask], y[mask]
    if xc.nunique() < 2 or yc.nunique() < 2:
        return (np.nan, np.nan, np.nan, np.nan)
    r, p = stats.pearsonr(xc, yc)
    rho, rp = stats.spearmanr(xc, yc)
    return (float(r), float(p), float(rho), float(rp))


def correlation_plot(
    matched_df: pd.DataFrame,
    score_columns: list[str],
    measured_col: str,
    out_png: Path,
    exon_col: Optional[str] = "Region (Exon)",
    include_combined: bool = True,
) -> dict:
    """Render per-exon regression panels with one line per score source.

    Scores are min-max normalized per panel to share a common x-axis;
    correlations are computed on raw values. When `include_combined=True`
    and multiple exons are present, an additional "all" panel pools the
    samples. Returns `{panel: {source: (pearson_r, pearson_p, spearman_rho, spearman_p)}}`.
    """
    import matplotlib.pyplot as plt
    import seaborn as sns

    if exon_col and exon_col in matched_df.columns:
        exons = sorted(matched_df[exon_col].dropna().unique().tolist())
        panels = [(str(e), matched_df[matched_df[exon_col] == e].copy()) for e in exons]
        if include_combined and len(exons) > 1:
            panels.append(("all", matched_df.copy()))
    else:
        panels = [("all", matched_df.copy())]
    if not panels:
        raise ValueError("No data to plot.")

    fig, axes = plt.subplots(1, len(panels), figsize=(7 * len(panels), 6), squeeze=False)
    palette = ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple"]
    stats_out: dict[str, dict[str, tuple[float, float, float, float]]] = {}

    for ax, (panel_name, sub) in zip(axes[0], panels):
        stats_out[panel_name] = {}
        if len(sub) < 3:
            ax.set_title(f"{panel_name} (n={len(sub)}, too few)")
            continue

        for col in score_columns:
            if col not in sub.columns:
                continue
            lo, hi = sub[col].min(), sub[col].max()
            sub[f"{col}_norm"] = (sub[col] - lo) / (hi - lo) if hi > lo else 0.0

        for i, col in enumerate(score_columns):
            if col not in sub.columns:
                continue
            r, p, rho, rp = _safe_corr(sub[col], sub[measured_col])
            stats_out[panel_name][col] = (r, p, rho, rp)
            label = (
                f"{col} (r={r:.3f}, ρ={rho:.3f})" if not np.isnan(r) else f"{col} (n/a)"
            )
            sns.regplot(
                data=sub,
                x=f"{col}_norm",
                y=measured_col,
                ax=ax,
                color=palette[i % len(palette)],
                label=label,
                scatter_kws={"alpha": 0.7, "edgecolors": "k", "linewidths": 0.5},
            )

        ax.set_xlabel("Predicted score (normalized)")
        ax.set_ylabel(measured_col)
        ax.set_title(f"Exon {panel_name} (n={len(sub)})")
        ax.legend(fontsize=8)

    plt.tight_layout()
    out_png = Path(out_png)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=150)
    plt.close(fig)
    return stats_out


def print_correlation_table(stats: dict) -> None:
    """Pretty-print the correlation results as a table."""
    def _fmt(x: float) -> str:
        return "n/a" if np.isnan(x) else f"{x:.4f}"

    rows = [
        {
            "Exon": exon,
            "Source": src,
            "Pearson r": _fmt(r),
            "Pearson p": _fmt(p),
            "Spearman ρ": _fmt(rho),
            "Spearman p": _fmt(rp),
        }
        for exon, by_src in stats.items()
        for src, (r, p, rho, rp) in by_src.items()
    ]
    if not rows:
        print("No correlation results.")
        return
    print("\n=== Correlation summary ===")
    print(pd.DataFrame(rows).to_string(index=False))
