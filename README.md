# AI4VPrimer

Automated degenerate primer design for viral surveillance.

AI4VPrimer designs PCR primers from aligned viral sequences, handling sequence diversity through IUPAC ambiguity codes. It selects optimal primer pairs for full-length sequencing and qPCR detection — all configurable through a no-code web interface.

---

## How It Works

```
                         ┌──────────────────────┐
                         │   Aligned FASTA file  │
                         └──────────┬───────────┘
                                    │
                                    ▼
                    ┌───────────────────────────────┐
                    │  ≤1000 seqs?       >1000 seqs?│
                    │  Standard mode     Subsample  │
                    │  (use all)         (1000 × N) │
                    └───────────────┬───────────────┘
                                    │
                                    ▼
                    ┌───────────────────────────────┐
                    │   Find Conserved Regions      │
                    │   (≥90% nucleotide identity)  │
                    └───────────────┬───────────────┘
                                    │
                                    ▼
                    ┌───────────────────────────────┐
                    │   Design Degenerate Primers   │
                    │   (IUPAC codes for diversity) │
                    └───────────────┬───────────────┘
                                    │
                                    ▼
                    ┌───────────────────────────────┐
                    │   Filter Primers              │
                    │   • Melting temperature (Tm)  │
                    │   • Primer dimer check        │
                    │   • GC content                │
                    │   • Degeneracy limit          │
                    └───────────────┬───────────────┘
                                    │
                                    ▼
                    ┌───────────────────────────────┐
                    │   Validate on Full Dataset    │
                    │   (if subsampled)             │
                    └───────────────┬───────────────┘
                                    │
                         ┌──────────┴──────────┐
                         ▼                     ▼
              ┌─────────────────┐   ┌─────────────────┐
              │  Sequencing     │   │  qPCR           │
              │  500-1200 bp    │   │  80-350 bp      │
              │  2 overlapping  │   │  1 best pair    │
              │  primer pairs   │   │                 │
              └─────────────────┘   └─────────────────┘
                         │                     │
                         └──────────┬──────────┘
                                    ▼
                         ┌─────────────────────┐
                         │  Output Report (.md) │
                         └─────────────────────┘
```

---

## Quick Start

```bash
git clone https://github.com/YOUR_USERNAME/AI4VPrimer.git
cd AI4VPrimer
pip install -r requirements.txt
python app.py
```

The browser opens automatically. Fill in your file paths, click Save, click Run.

---

## Files

| File | Purpose |
|------|---------|
| `app.py` | Web server and API |
| `index.html` | Browser interface |
| `primer_pipeline.py` | Core pipeline |
| `requirements.txt` | Python dependencies |
| `README.md` | This file |
| `GUIDE.md` | Detailed usage guide |

---

## Requirements

- Python 3.7+
- Flask ≥ 3.0
- Biopython, NumPy, primer3-py

```bash
pip install -r requirements.txt
```

---

## Input

An **aligned** FASTA file (use MAFFT or MUSCLE to align beforehand).

## Output

A markdown report with:
- All candidate primers (Tm, GC%, sensitivity)
- Recommended sequencing primer pairs (overlapping amplicons)
- Recommended qPCR primer pair
- Reverse complements for PCR ordering

