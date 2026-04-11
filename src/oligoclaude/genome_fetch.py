"""GRCh38 FASTA cache — download once, reuse across configs.

Streams the gzipped GENCODE primary assembly through gunzip directly into
the cache file so we never keep the ~900 MB .gz and the ~3 GB .fa on disk
simultaneously. The final uncompressed file is unavoidable because pyfaidx
needs random access.
"""
from __future__ import annotations

import gzip
import shutil
import sys
import urllib.request
from pathlib import Path
from typing import Optional

CACHE_DIR = Path.home() / ".oligoclaude" / "genomes"
HG38_FILENAME = "GRCh38.primary_assembly.genome.fa"
HG38_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/"
    "release_46/GRCh38.primary_assembly.genome.fa.gz"
)


def default_hg38_path() -> Path:
    return CACHE_DIR / HG38_FILENAME


def _iter_with_progress(src, total: Optional[int], label: str):
    """Copy src to a destination, printing percentage progress."""
    chunk = 1024 * 1024
    written = 0
    last_pct = -1
    while True:
        buf = src.read(chunk)
        if not buf:
            break
        written += len(buf)
        if total:
            pct = int(written * 100 / total)
            if pct != last_pct:
                sys.stderr.write(f"\r  {label}: {pct}%")
                sys.stderr.flush()
                last_pct = pct
        yield buf
    if total:
        sys.stderr.write("\n")


def ensure_hg38_fasta(
    cache_dir: Optional[Path] = None, *, verbose: bool = True
) -> Path:
    """Download and stream-decompress GRCh38 if not already cached.

    Returns the path to the uncompressed FASTA. Idempotent: returns the
    existing file immediately if it looks usable (> 2 GB).
    """
    cache_dir = Path(cache_dir) if cache_dir else CACHE_DIR
    cache_dir.mkdir(parents=True, exist_ok=True)
    dst = cache_dir / HG38_FILENAME

    if dst.exists() and dst.stat().st_size > 2 * 1024**3:
        if verbose:
            print(f"Using cached GRCh38 FASTA: {dst}")
        return dst

    tmp = dst.with_suffix(dst.suffix + ".partial")
    if verbose:
        print(f"Downloading GRCh38 primary assembly FASTA...")
        print(f"  Source: {HG38_URL}")
        print(f"  Destination: {dst}")
        print("  Streaming gunzip (~900 MB download, ~3 GB on disk)...")

    req = urllib.request.Request(HG38_URL, headers={"User-Agent": "oligoclaude"})
    with urllib.request.urlopen(req) as resp:
        total = resp.headers.get("Content-Length")
        total_int = int(total) if total and total.isdigit() else None
        with gzip.GzipFile(fileobj=resp) as gz, open(tmp, "wb") as out:
            if verbose and total_int:
                for chunk in _iter_with_progress(gz, total_int, "download (gz)"):
                    out.write(chunk)
            else:
                shutil.copyfileobj(gz, out, length=1024 * 1024)

    tmp.replace(dst)
    if verbose:
        size_gb = dst.stat().st_size / 1024**3
        print(f"  Decompressed FASTA: {size_gb:.2f} GB → {dst}")
    return dst


def resolve_fasta_path(
    cfg_value: Optional[Path], *, auto_fetch: bool = True, verbose: bool = True
) -> Path:
    """Resolve the FASTA path, auto-downloading the default cache if needed.

    - If `cfg_value` is given and exists, return it.
    - If `cfg_value` is given but missing, fall through to the cache.
    - If nothing exists and `auto_fetch` is True, download to the cache.
    """
    if cfg_value and Path(cfg_value).exists():
        return Path(cfg_value)

    cached = default_hg38_path()
    if cached.exists() and cached.stat().st_size > 2 * 1024**3:
        if verbose and cfg_value:
            print(
                f"fasta_path {cfg_value} not found — falling back to "
                f"cached {cached}"
            )
        return cached

    if auto_fetch:
        return ensure_hg38_fasta(verbose=verbose)

    raise FileNotFoundError(
        f"FASTA not found: {cfg_value}. Run `oligoclaude fetch-genome` "
        "or download GRCh38.primary_assembly.genome.fa manually."
    )
