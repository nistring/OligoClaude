"""Sequence utilities, ASO enumeration, and experimental-data helpers."""
from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
_BASE_IDX = {"A": 0, "C": 1, "G": 2, "T": 3}
_MATCH_SUFFIX_RE = re.compile(r"_m\d+$")


def reverse_complement(seq: str) -> str:
    return seq.translate(COMP)[::-1]


def load_reference_sequence(fasta_path: Path, chrom: str, start: int, end: int) -> str:
    """Load a genomic subsequence (0-based, end-exclusive) via pyfaidx."""
    from pyfaidx import Fasta

    fa = Fasta(str(fasta_path), as_raw=True, sequence_always_upper=True)
    if chrom not in fa:
        raise KeyError(
            f"Chromosome {chrom!r} not found in {fasta_path}. "
            f"Available: {list(fa.keys())[:5]}..."
        )
    return str(fa[chrom][start:end])


def one_hot_encode(seq: str) -> np.ndarray:
    """One-hot encode a DNA sequence. N → [0.25, 0.25, 0.25, 0.25]."""
    arr = np.zeros((len(seq), 4), dtype=np.float32)
    for i, b in enumerate(seq.upper()):
        idx = _BASE_IDX.get(b)
        if idx is None:
            arr[i] = 0.25
        else:
            arr[i, idx] = 1.0
    return arr


@dataclass
class AsoCandidate:
    """A single ASO candidate.

    `position` is the 0-based offset within `ref_seq` (the reference region
    loaded by the workflow); downstream code anchors BED/scoring windows at
    ref_seq[0]. `length` may vary in experimental mode.
    """

    aso_id: str
    aso_sequence_antisense: str
    genomic_target_seq: str
    position: int
    length: int
    measured: Optional[float] = None
    exon_label: Optional[str] = None


def _strip_match_suffix(aso_id: str) -> str:
    return _MATCH_SUFFIX_RE.sub("", aso_id)


def enumerate_sliding(
    ref_seq: str, start_rel: int, end_rel: int, aso_length: int, step: int = 1
) -> list[AsoCandidate]:
    """Sliding-window ASO enumeration across [start_rel, end_rel) in ref_seq."""
    cands: list[AsoCandidate] = []
    for i in range(start_rel, end_rel - aso_length + 1, step):
        target = ref_seq[i : i + aso_length]
        cands.append(
            AsoCandidate(
                aso_id=f"win_{i - start_rel}",
                aso_sequence_antisense=reverse_complement(target),
                genomic_target_seq=target,
                position=i,
                length=aso_length,
            )
        )
    return cands


def enumerate_from_experimental(
    exp_df: pd.DataFrame,
    ref_seq: str,
    variant_interval_start_rel: int,
    variant_interval_end_rel: int,
    aso_seq_col: str = "ASO sequence",
    aso_id_col: str = "ASO_ID",
    measured_col: Optional[str] = "Measured (RT-PCR)",
    exon_col: Optional[str] = "Region (Exon)",
) -> list[AsoCandidate]:
    """Build candidates from an experimental CSV by reverse-complement lookup.

    Each experimental ASO sequence (antisense) is reverse-complemented and
    searched within [variant_interval_start_rel, variant_interval_end_rel)
    in ref_seq. Every occurrence yields one AsoCandidate; multiple matches
    get `_m2`, `_m3`, … suffixes. `position` is the ref_seq offset.
    """
    cands: list[AsoCandidate] = []
    for idx, row in exp_df.iterrows():
        antisense = str(row[aso_seq_col]).strip().upper()
        if not antisense:
            continue
        target = reverse_complement(antisense)
        length = len(target)
        exp_id = str(row[aso_id_col]) if aso_id_col in row else f"exp_{idx}"
        measured = (
            float(row[measured_col])
            if measured_col and measured_col in row and pd.notna(row[measured_col])
            else None
        )
        exon_label = (
            str(row[exon_col])
            if exon_col and exon_col in row and pd.notna(row[exon_col])
            else None
        )

        search_from = variant_interval_start_rel
        match_idx = 0
        while True:
            pos = ref_seq.find(target, search_from, variant_interval_end_rel)
            if pos == -1 or pos + length > variant_interval_end_rel:
                break
            match_idx += 1
            cands.append(
                AsoCandidate(
                    aso_id=exp_id if match_idx == 1 else f"{exp_id}_m{match_idx}",
                    aso_sequence_antisense=antisense,
                    genomic_target_seq=target,
                    position=pos,
                    length=length,
                    measured=measured,
                    exon_label=exon_label,
                )
            )
            search_from = pos + 1
    return cands


def load_experimental(csv_path: Path) -> pd.DataFrame:
    """Load and whitespace-normalize an experimental CSV."""
    df = pd.read_csv(csv_path)
    df.columns = df.columns.str.strip()
    return df


def match_scores_to_experimental(
    exp_df: pd.DataFrame,
    candidates: list[AsoCandidate],
    scores: dict[str, np.ndarray],
    pred_length: Optional[int] = None,
) -> pd.DataFrame:
    """Attach per-source predicted scores to each experimental row.

    For each row, reverse-complement its ASO sequence, slide a window of
    `pred_length` across the target, look up each k-mer in the candidate
    score table, and average matches. Unmatched rows get NaN.
    """
    if pred_length is None and candidates:
        pred_length = candidates[0].length

    lookups = {
        src: {c.genomic_target_seq: float(val) for c, val in zip(candidates, arr)}
        for src, arr in scores.items()
    }

    out = exp_df.copy()
    for src in scores:
        out[src] = np.nan

    for idx, row in out.iterrows():
        antisense = str(row.get("ASO sequence", "")).strip().upper()
        if not antisense:
            continue
        target = reverse_complement(antisense)
        windows = [target[i : i + pred_length] for i in range(len(target) - pred_length + 1)]
        for src, lookup in lookups.items():
            vals = [lookup[w] for w in windows if w in lookup]
            if vals:
                out.at[idx, src] = float(np.mean(vals))
    return out


def aggregate_experimental_candidates(
    candidates: list[AsoCandidate],
    scores: dict[str, np.ndarray],
    measured_col: str = "Measured (RT-PCR)",
    aso_id_col: str = "ASO_ID",
    exon_col: str = "Region (Exon)",
) -> pd.DataFrame:
    """Collapse multi-match experimental candidates to one row per base ASO.

    Candidates with IDs like `a1` and `a1_m2` (same experimental ASO, two
    hits in ref_seq) are grouped by the stripped base ID and scores are
    averaged. Returns [aso_id_col, exon_col, measured_col, "ASO sequence",
    *score_keys].
    """
    if not candidates:
        return pd.DataFrame()

    rows = [
        {
            "_base": _strip_match_suffix(c.aso_id),
            aso_id_col: _strip_match_suffix(c.aso_id),
            exon_col: c.exon_label,
            measured_col: c.measured,
            "ASO sequence": c.aso_sequence_antisense,
            **{k: float(scores[k][i]) for k in scores},
        }
        for i, c in enumerate(candidates)
    ]
    df = pd.DataFrame(rows)
    agg = {aso_id_col: "first", exon_col: "first", measured_col: "first", "ASO sequence": "first"}
    agg.update({k: "mean" for k in scores})
    return df.groupby("_base", as_index=False, sort=False).agg(agg).drop(columns="_base")
