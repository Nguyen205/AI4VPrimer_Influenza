# Degenerate Primer Design Pipeline

Automated pipeline for designing degenerate PCR primers from sequence alignments, with selection for full-length sequencing and qPCR applications.

## Overview

This pipeline:
1. Identifies conserved regions from aligned sequences
2. Designs degenerate primers with IUPAC ambiguity codes
3. Filters by Tm, homodimer, GC content, and consecutive bases
4. Selects up to 3 primers per conserved region
5. Recommends optimal primer pairs for sequencing (500-1200 bp) and qPCR (80-350 bp)
6. Validates sensitivity on cross-panels (optional)
7. Calculates specificity against non-target sequences (optional)

## Requirements

```bash
pip install biopython numpy primer3-py
```

## Files

| File | Purpose |
|------|---------|
| `primer_pipeline.py` | Main design & selection pipeline |
| `sensitivity and specificity.ipynb` | Additional validation of primers |
| `merge fasta.ipynb` | Combine multiple FASTA files |

---

## Usage

### Step 1: Configure the Pipeline

Edit the CONFIG section in `primer_pipeline.py`:

```python
GENE_NAME  = 'H5_VN'
ALN_FILE   = '/path/to/aligned_sequences.fasta'  # Required, must be aligned

# Cross-panels for sensitivity check (optional, set to None or {} if not needed)
CROSS_PANELS = {
    'Asia':   '/path/to/asia_panel.fasta',
    'Global': '/path/to/global_panel.fasta',
}

# Specificity check against non-target sequences (optional, set to None if not needed)
SPECIFICITY_FILE = '/path/to/non_target_sequences.fasta'

OUT_FILE = '/path/to/output.md'
PRIMERS_PER_REGION = 3  # Number of primers to keep per conserved region
```

### Step 2: Run the Pipeline

```bash
python primer_pipeline.py
```

### Step 3: Additional Validation (Optional)

Use `sensitivity and specificity.ipynb` for further primer evaluation:

```python
calculate_metrics(fasta_file1, FP, RP, merged_fasta_file)
```

---

## Pipeline Workflow

```
Aligned FASTA → Conserved Regions → Primer Construction → Trimming
                                                              ↓
                              Filtering ← Degeneracy Reduction ← End Fixing
                                  ↓
                    Tiebreaking (GC content, no AAAA/CCCC/GGGG/TTTT)
                                  ↓
                    Keep up to 3 primers per region
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
              Sensitivity (cross-panels) & Specificity (non-targets)
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
| `PRIMERS_PER_REGION` | 3 | Max primers kept per conserved region |

---

## Tiebreaker Logic

When multiple primers have the same sensitivity, the pipeline selects based on:
1. **GC content** closer to 50% (ideal: 40-60%)
2. **No consecutive bases** (avoids AAAA, CCCC, GGGG, TTTT)

---

## Sensitivity vs Specificity

- **Sensitivity:** % of target sequences where BOTH primers match the same sequence
- **Specificity:** % of non-target sequences where the primer combo does NOT match

---

## Target Selection & Sequence Data

After inspecting nucleotide identity across viral genomes, we focused on **M**, **NA**, and **HA** segments for degenerate primer design.

Sequences in this repository include:
- Segment sequences (M, NA, HA) from various **Influenza A subtypes** (H1N2, H7N9, etc.)
- Genomes from other **aerosolized viruses** (PRCV, Norovirus)
- Merged FASTA files for sensitivity/specificity calculations

---

## Output

- **Primer Pool:** All primers passing filters with location, Tm, GC%, sensitivity
- **Sequencing Selection:** Two overlapping primer pairs with true PCR sensitivity and specificity
- **qPCR Selection:** Single best primer pair with sensitivity and specificity

**Note:** Reverse primers are in FASTA orientation. Use the reverse complement (RC) for PCR ordering.

---

## Notes

- **Input alignment:** Must be pre-aligned (use MAFFT, MUSCLE, etc.)
- **Cross-panels:** Optional, can be unaligned raw FASTA
- **Specificity file:** Optional, contains non-target sequences
- **Degeneracy:** Fwd × Rev product should be ≤24 (preferred) or ≤32 (max)
