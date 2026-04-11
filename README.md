# OligoClaude

Predict antisense oligonucleotide (ASO) efficacy using **AlphaGenome** and **SpliceAI**, compare to experimental RT-PCR data, and export UCSC-ready BED tracks — all from a single JSON config.

Exposed three ways:
1. CLI (`oligoclaude run --config …`)
2. Python library (`from oligoclaude import run_workflow`)
3. **Local MCP server** so Claude Desktop / Claude Code can call it as a connector

## Install

### Linux / macOS

```bash
git clone <this-repo>
cd OligoClaude
pip install -e ".[spliceai]"          # AlphaGenome + SpliceAI
pip install -e .                       # AlphaGenome only (skip SpliceAI)
```

### Windows

The SpliceAI stack (`openspliceai`) transitively depends on `pysam`, which has
no Windows wheels. Install the SpliceAI model code without its pysam-dependent
CLI:

```powershell
git clone <this-repo>
cd OligoClaude
pip install openspliceai --no-deps     # pulls the pure-torch model class only
pip install -e .                       # all other deps
```

OligoClaude only imports the pysam-free path
`openspliceai.train_base.openspliceai.SpliceAI`, so `--no-deps` is enough.

For unit tests: `pip install -e .[dev]`

### One-time setup

```bash
oligoclaude set-api-key YOUR_ALPHAGENOME_KEY   # saved to ~/.oligoclaude/credentials.json (0600)
oligoclaude fetch-genome                        # downloads GRCh38 FASTA to ~/.oligoclaude/genomes/ (~3 GB)
oligoclaude fetch-spliceai-weights              # downloads MANE-10000nt ensemble (~40 MB × 5)
```

All three commands are idempotent. The genome and SpliceAI weights are
cached under `~/.oligoclaude/` and reused across configs. The genome fetcher
streams the gzipped download directly through gunzip so you never store both
the `.gz` and the uncompressed `.fa` simultaneously.

You can also set the API key via the `ALPHAGENOME_API_KEY` environment
variable instead of the credentials file.

## Config

Create a JSON config (see `config/SETD5_e1.json` for a full example):

```json
{
  "gene_symbol": "SETD5",
  "assembly": "hg38",
  "gtf_url": "https://storage.googleapis.com/alphagenome/reference/gencode/hg38/gencode.v46.annotation.gtf.gz.feather",
  "results_dir": "results",
  "data_dir": "data",
  "track_filter": "polyA plus RNA-seq",
  "ontology_terms": ["CL:0000127", "CL:0002319"],
  "requested_outputs": ["RNA_SEQ", "SPLICE_SITE_USAGE"],
  "exon_intervals": [9429826, 9430051],
  "flank": [200, 200],
  "strand": "+",
  "ASO_length": 18,
  "aso_step": 1,
  "experimental_data": "ASOseq_incl.csv",
  "target_mode": "exclude",
  "spliceai_batch": 12
}
```

Required fields: `gene_symbol`, `assembly`, `gtf_url`, `exon_intervals`, `strand`, `ASO_length`, `flank`, `ontology_terms`, `requested_outputs`.

**Do NOT put `dna_api_key` in the config file.** Use `oligoclaude set-api-key`
or the `ALPHAGENOME_API_KEY` environment variable. A legacy `dna_api_key`
field is still honored for backward compatibility, but triggers a
`DeprecationWarning` because it gets committed to version control.

**`fasta_path` is optional.** If omitted, OligoClaude uses the auto-downloaded
GRCh38 FASTA at `~/.oligoclaude/genomes/GRCh38.primary_assembly.genome.fa`
(run `oligoclaude fetch-genome` once, or let the first `run` invocation
fetch it automatically).

Optional fields:
- `experimental_data`: path to a CSV with columns `ASO_ID`, `ASO sequence`, `Measured (RT-PCR)`, optionally `Region (Exon)`. If present, OligoClaude scores those ASOs directly and produces a correlation plot; otherwise it runs a sliding window across `[exon_start - flank[0], exon_end + flank[1]]`.
- `aso_step`: sliding window step size (default 1).
- `target_mode`: `"exclude"` (default) or `"include"`. Only affects interpretation, not raw scores.
- `spliceai_batch`, `spliceai_threads`: SpliceAI CPU batching knobs.

Relative paths inside the JSON are resolved against the config file's directory.

## CLI usage

```bash
# Full workflow (AlphaGenome + SpliceAI + correlation plot)
oligoclaude run --config config/SETD5_e1.json -v

# Skip SpliceAI (faster; AlphaGenome only)
oligoclaude run --config config/SETD5_e1.json --skip-spliceai

# Skip AlphaGenome (no API key needed)
oligoclaude run --config config/SETD5_e1.json --skip-alphagenome
```

Outputs land in `<results_dir>/ASO/`:
- `<config_name>_ASO_scores.csv` — raw scores per ASO, all sources
- `<config_name>_ASO_<source>.bed`, `<config_name>_ASO_<source>_full.bed` — UCSC custom track files
- `<config_name>_correlation.png` — if experimental data provided
- `<config_name>_experimental_matched.csv` — aggregated predictions vs measured

## UCSC custom track upload

OligoClaude writes BED files matching the format used by `genome.ucsc.edu`'s custom tracks. To view them:

1. Go to `https://genome.ucsc.edu/cgi-bin/hgCustom?db=hg38`
2. Click **Choose File** and pick one of the generated `*.bed` files
3. Click **Submit**

The tool prints these instructions (with the exact file paths) at the end of every run.

## Use from Claude (MCP connector)

OligoClaude ships a local **MCP server** that Claude Code and Claude Desktop can call as a first-class connector, exposing a single tool: **`predict_aso_efficacy`**.

### Register in Claude Code

```bash
pip install -e /path/to/OligoClaude
claude mcp add oligoclaude -- oligoclaude-mcp
```

### Register in Claude Desktop

Edit `~/.config/Claude/claude_desktop_config.json` (Linux),
`~/Library/Application Support/Claude/claude_desktop_config.json` (macOS), or
`%APPDATA%\Claude\claude_desktop_config.json` (Windows):

```json
{
  "mcpServers": {
    "oligoclaude": {
      "command": "oligoclaude-mcp",
      "env": {
        "ALPHAGENOME_API_KEY": "your-alphagenome-key-here"
      }
    }
  }
}
```

Alternatively, omit the `env` block and run `oligoclaude set-api-key` once —
the server reads `~/.oligoclaude/credentials.json` on startup.

Restart Claude Desktop. Then in any conversation you can say:

> *"Use the oligoclaude MCP to predict ASO efficacy for /home/me/OligoClaude/config/SETD5_e1.json"*

Claude will invoke `predict_aso_efficacy(config_path=...)` and receive a dict with the scores CSV path, BED file paths, correlation plot path, per-exon Pearson/Spearman stats, and UCSC upload instructions.

### Tool signature

```python
predict_aso_efficacy(
    config_path: str,
    skip_alphagenome: bool = False,
    skip_spliceai: bool = False,
    samples_max: int = 20,
) -> dict
```

## How scoring works

**AlphaGenome** — for each ASO position, the target nucleotides are replaced with `N`s and the resulting sequence is predicted via `model.predict_sequences`. Per-output-type efficacy follows the `diff_mean` formula from AlphaGenome_ASO/aso.ipynb:

```
score = alt[exon].mean() / ref[gene_body].mean() * alt[gene_body].mean() - ref[exon].mean()
```

**SpliceAI** — the reference sequence (`SL=5000` centered on the exon, plus `CL_max=10000` flanking context) is one-hot encoded once; per-ASO variants are built by `clone() + slice-assign` to set `[0.25, 0.25, 0.25, 0.25]` at the ASO position. The 5-model OpenSpliceAI MANE-10000nt ensemble runs on CPU with configurable batching. The efficacy score uses the same `diff_mean` formula applied to per-position donor+acceptor probability sums.

Both scores are positive when the ASO weakens exon usage (exon skipping) and negative when it strengthens it.

## Library usage

```python
from oligoclaude import run_workflow

result = run_workflow("config/SETD5_e1.json", verbose=True)
print(result.scores_csv)
print(result.stats)
```

## License

See `LICENSE`.
