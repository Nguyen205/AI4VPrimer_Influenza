# Degenerate Primer Design Pipeline

Automated pipeline for designing degenerate PCR primers from sequence alignments, with selection for full-length sequencing and qPCR applications.

## Overview

This pipeline:
1. Identifies conserved regions from aligned sequences
2. Designs degenerate primers with IUPAC ambiguity codes
3. Filters by Tm and homodimer formation
4. Selects optimal primer pairs for sequencing (500-1200 bp) and qPCR (80-350 bp)
5. Validates sensitivity/specificity against target and non-target sequences

## Requirements

```bash
pip install biopython numpy primer3-py
```

## Files

| File | Purpose |
|------|---------|
| `primer_pipeline.py` | Main design & selection pipeline |
| `sensitivity and specificity.ipynb` | Validate primers against target/non-target sequences |
| `merge fasta.ipynb` | Combine multiple FASTA files for analysis |

---

## Usage

### Step 1: Design Primers

Edit the CONFIG section in `primer_pipeline.py`:

```python
GENE_NAME  = 'H5_VN'
ALN_FILE   = '/path/to/aligned_sequences.fasta'  # Must be aligned
CROSS_PANELS = {
    'Asia':   '/path/to/asia_panel.fasta',       # Can be unaligned
    'Global': '/path/to/global_panel.fasta',
}
OUT_FILE   = '/path/to/output.md'
```

Run:
```bash
python primer_pipeline.py
```

### Step 2: Validate Sensitivity & Specificity

Use `sensitivity and specificity.ipynb` to evaluate primer performance:

```python
calculate_metrics(fasta_file1, FP, RP, merged_fasta_file)
```

- **Sensitivity:** % of target sequences detected by both primers
- **Specificity:** % of non-target sequences correctly excluded

### Step 3: Merge FASTA Files (if needed)

Use `merge fasta.ipynb` to consolidate sequence files:

```python
merge_fasta_files(file_list, output_file)
```

---

## Pipeline Workflow

```
Aligned FASTA → Conserved Regions → Primer Construction → Trimming
                                                              ↓
                              Filtering ← Degeneracy Reduction ← End Fixing
                                  ↓
                           Primer Pool
                                  ↓
                    ┌─────────────┴─────────────┐
                    ↓                           ↓
            Sequencing Selection          qPCR Selection
            (500-1200 bp, overlap)        (80-350 bp)
                    ↓                           ↓
                    └─────────────┬─────────────┘
                                  ↓
                    Sensitivity/Specificity Validation
```

---

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `IDENTITY_THRESHOLD` | 90% | Min conservation for region selection |
| `MIN_REGION_LEN` | 16 bp | Min conserved region length |
| `MAX_PRIMER_LEN` | 24 bp | Max primer length |
| `MIN_TM` | 51.4°C | Min melting temperature |
| `SENS_TARGET` | 80% | Target true PCR sensitivity |
| `COMBO_DEGEN_MAX` | 32 | Max combined degeneracy (fwd × rev) |

---

## Target Selection & Sequence Data

After inspecting nucleotide identity across viral genomes, we focused on **M**, **NA**, and **HA** segments for degenerate primer design.

Sequences in this repository include:
- Segment sequences (M, NA, HA) from various **Influenza A subtypes** (H1N2, H7N9, etc.)
- Genomes from other
 **aerosolized viruses** (PRCV, Norovirus)
- Merged FASTA files for sensitivity/specificity calculations

---

## Output

### Primer Pool
All primers passing filters with location, sequence, Tm, and sensitivity on each panel.

### Selected Primers
- **Sequencing:** Two overlapping primer pairs (4 primers total)
- **qPCR:** Single best primer pair

**Note:** Reverse primers are shown in FASTA orientation. Use the reverse complement (RC) for PCR ordering.

---

## Notes

- **Input alignment:** Must be pre-aligned (use MAFFT, MUSCLE, etc.)
- **Cross-panels:** Can be unaligned raw FASTA
- **True PCR sensitivity:** Both forward AND reverse must match the same sequence
- **Degeneracy:** Fwd × Rev product should be ≤24 (preferred) or ≤32 (max)


