# OligoMCP

Design and score antisense oligonucleotides (ASOs) with **AlphaGenome**
and **SpliceAI**. Drive it from natural language in Claude (MCP) or
from the CLI.

## Install

```bash
pip install -e .
oligomcp init        # prompts for AlphaGenome key
```

SpliceAI weights ship with the package. Reference sequence is fetched
from UCSC on demand — no genome download required.

## Use from Claude

Register once:

```bash
claude mcp add oligomcp -- oligomcp-mcp
```

Then just ask:

> *"Design ASOs to skip SETD5 exon 11 and show the top 10 by SpliceAI."*
>
> *"Score ASOs for SMN2 exon 7 in motor-neuron tracks."*
>
> *"Run `config/SETD5_e1.json` and correlate with the experimental CSV."*

Claude picks the right tool and returns scores, BED files, and plots.

## Use from the CLI

```bash
oligomcp run --config config/SETD5_e1.json -v
oligomcp run --config config/SETD5_e1.json --skip-alphagenome    # SpliceAI only
```

Minimal config:

```json
{ "gene_symbol": "SETD5", "exon_intervals": [9429826, 9430051] }
```

See `config/SETD5_e1.json` for a full example.

## Designing against a patient genome (variants)

Add a `variants` list to edit mutations into the reference before ASO
enumeration — useful when a known SNP would disrupt ASO affinity or
when the ASO needs to overlap a disease-causing variant. SNVs **and
indels** are supported. Coordinates are always reference-genome (as in
VCF).

```json
{
  "gene_symbol": "SMN2",
  "exon_intervals": [70070640, 70070751],
  "variants": [
    "chr5:70070740:G:A",
    "rs1800112",
    "c.840C>T",
    { "id": "my_del", "chrom": "chr5", "position": 70070700,
      "ref": "AAG", "alt": "" }
  ]
}
```

Accepted formats (mix freely): VCF-style, HGVS `g.`, HGVS `c.`, rsID,
ClinVar `VCV…`, or an explicit dict.

## Outputs

Land in `results/<gene_symbol>/<timestamp>/` — one folder per run.

| File | What |
|---|---|
| `<name>_scores.csv` | Per-ASO scores |
| `<name>_<source>.bed` | UCSC/IGV track (red = inclusion, blue = skipping) |
| `<name>_correlation.png` | Predicted ↔ measured regression |

Positive score = strengthens inclusion. Negative = promotes skipping.
For `target_mode="exclude"`, pick the most negative.

## License

See `LICENSE`.
