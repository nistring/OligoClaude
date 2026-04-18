"""Tests for BED export."""
from pathlib import Path

import numpy as np

from oligomcp.core import AsoCandidate
from oligomcp.output import write_bed


def _mk_candidates(n: int, length: int = 18) -> list:
    cands = []
    for i in range(n):
        cands.append(
            AsoCandidate(
                aso_id=f"win_{i}",
                aso_sequence_antisense="N" * length,
                genomic_target_seq="N" * length,
                position=i,
                length=length,
            )
        )
    return cands


def _write(tmp_path: Path, cands, scores, **kwargs) -> Path:
    defaults = dict(
        results_dir=tmp_path,
        config_name="TEST",
        source="RNA_SEQ",
        chrom="chr1",
        strand="+",
        candidates=cands,
        scores=scores,
        variant_interval_start=1000,
    )
    defaults.update(kwargs)
    return write_bed(**defaults)


def test_write_bed_header_and_columns(tmp_path: Path):
    cands = _mk_candidates(5)
    scores = np.array([0.5, -0.3, 0.1, -0.2, 0.4], dtype=np.float32)

    path = _write(tmp_path, cands, scores)
    assert path.name == "TEST_RNA_SEQ.bed"

    lines = path.read_text().strip().split("\n")
    assert lines[0].startswith("track name=TEST_RNA_SEQ ")
    assert "useScore=1" in lines[0]
    assert 'visibility="pack"' in lines[0]

    # all 5 candidates (3 positive + 2 negative)
    assert len(lines[1:]) == 5


def test_write_bed_strand_inversion_plus(tmp_path: Path):
    cands = _mk_candidates(1)
    scores = np.array([0.5], dtype=np.float32)
    path = _write(tmp_path, cands, scores, source="X", config_name="T", variant_interval_start=0)
    cols = path.read_text().strip().split("\n")[1].split("\t")
    assert cols[5] == "-"


def test_write_bed_strand_inversion_minus(tmp_path: Path):
    cands = _mk_candidates(1)
    scores = np.array([0.5], dtype=np.float32)
    path = _write(tmp_path, cands, scores, source="X", config_name="T",
                  strand="-", variant_interval_start=0)
    cols = path.read_text().strip().split("\n")[1].split("\t")
    assert cols[5] == "+"


def test_write_bed_color_positive_is_red(tmp_path: Path):
    cands = _mk_candidates(1)
    scores = np.array([0.5], dtype=np.float32)
    path = _write(tmp_path, cands, scores, source="X", config_name="T", variant_interval_start=0)
    rgb = path.read_text().strip().split("\n")[1].split("\t")[8].split(",")
    assert rgb[0] == "255"


def test_write_bed_color_negative_is_blue(tmp_path: Path):
    cands = _mk_candidates(1)
    scores = np.array([-0.5], dtype=np.float32)
    path = _write(tmp_path, cands, scores, source="X", config_name="T", variant_interval_start=0)
    rgb = path.read_text().strip().split("\n")[1].split("\t")[8].split(",")
    assert rgb[2] == "255"


def test_write_bed_absolute_coordinates(tmp_path: Path):
    cands = _mk_candidates(1, length=18)
    scores = np.array([0.3], dtype=np.float32)
    path = _write(tmp_path, cands, scores, source="X", config_name="T",
                  chrom="chr3", variant_interval_start=9429626)
    cols = path.read_text().strip().split("\n")[1].split("\t")
    assert cols[0] == "chr3"
    assert int(cols[1]) == 9429626
    assert int(cols[2]) == 9429644
    assert int(cols[6]) == 9429626
    assert int(cols[7]) == 9429644


def test_write_bed_score_is_integer_0_to_1000(tmp_path: Path):
    """UCSC BED spec: score column must be an integer in [0, 1000]."""
    cands = _mk_candidates(5)
    scores = np.array([0.0005, -0.0006, 0.0003, -0.0002, 0.0004], dtype=np.float32)

    path = _write(tmp_path, cands, scores, source="SpliceAI", config_name="T",
                  variant_interval_start=0)
    data_lines = path.read_text().strip().split("\n")[1:]
    assert data_lines, "expected BED data rows"
    for line in data_lines:
        score_col = line.split("\t")[4]
        assert score_col.lstrip("-").isdigit(), f"score column is not integer: {score_col!r}"
        val = int(score_col)
        assert 0 <= val <= 1000, f"score {val} out of BED spec range 0..1000"


def test_write_bed_header_has_itemRgb_on(tmp_path: Path):
    """Without itemRgb="On" in the track header, UCSC ignores the per-row RGB."""
    cands = _mk_candidates(2)
    scores = np.array([0.5, -0.3], dtype=np.float32)
    path = _write(tmp_path, cands, scores, source="SpliceAI", config_name="T",
                  variant_interval_start=0)
    header = path.read_text().splitlines()[0]
    assert 'itemRgb="On"' in header, f"missing itemRgb=\"On\" in {path.name}: {header!r}"
