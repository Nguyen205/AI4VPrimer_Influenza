# Degenerate Primer Design & Selection Pipeline

A unified Python pipeline for designing degenerate PCR primers from multiple sequence alignments, with automatic selection of optimal primer pairs for full-length sequencing and qPCR applications.

## Features

- **Automated conserved region identification** from sequence alignments
- **Degenerate primer design** using IUPAC ambiguity codes
- **Thermodynamic filtering** (Tm, homodimer check)
- **Cross-panel sensitivity validation** against multiple sequence databases
- **Dual application support:**
  - Full-length sequencing (overlapping amplicons, 500-1200 bp)
  - qPCR (short amplicons, 80-350 bp)
- **Degeneracy control** at both individual primer and combo levels

## Requirements

### Dependencies

```bash
pip install biopython numpy primer3-py
```

### Input Files

1. **Aligned FASTA file** (required) - Multiple sequence alignment of target gene
2. **Cross-panel FASTA files** (optional) - Unaligned sequences for sensitivity validation

## Quick Start

### 1. Configure the Pipeline

Edit the `CONFIG` section at the top of `primer_pipeline.py`:

```python
# ══════════════════════════════════════════════════════════════════════════════
# CONFIG - Modify this section for your gene
# ══════════════════════════════════════════════════════════════════════════════
GENE_NAME  = 'H5_VN'
ALN_FILE   = '/path/to/your/aligned_sequences.fasta'
CROSS_PANELS = {
    'Asia':   '/path/to/asia_panel.fasta',
    'Global': '/path/to/global_panel.fasta',
}
OUT_FILE   = '/path/to/output/primer_output.md'
```

### 2. Run the Pipeline

```bash
python primer_pipeline.py
```

### 3. View Results

The pipeline outputs:
- Console progress and summary
- Markdown report with primer pool and recommendations

## How It Works

### Stage 1: Conserved Region Identification

Scans the alignment to find regions where nucleotide identity ≥90% (configurable). Allows bridging of 1 low-identity position to capture longer regions.

### Stage 2: Primer Construction

For each conserved region, builds a degenerate consensus primer:
- Positions with ≥90% identity → use dominant base
- Positions with minor variants ≥5% → use IUPAC ambiguity code

### Stage 3: Primer Optimization

1. **Trimming:** Finds optimal 16-24 bp window with best sensitivity
2. **Degeneracy reduction:** Removes unnecessary ambiguous bases (if sensitivity loss <2%)
3. **End fixing:** Ensures first/last 2 positions are unambiguous (two strategies: replace or extend)

### Stage 4: Filtering

Primers must pass:
- Tm ≥51.4°C (primer and reverse complement)
- Homodimer Tm ≤35°C and ΔG ≥ -9 kcal/mol

### Stage 5: Selection

**Full-length sequencing:**
- Selects two overlapping primer pairs (4 primers)
- Amplicon size: 500-1200 bp each
- Overlap: ≤200 bp (relaxed if needed)
- Target: ≥80% true PCR sensitivity

**qPCR:**
- Selects single best primer pair
- Amplicon size: 80-350 bp
- Ranked by true PCR sensitivity

## Parameters

All parameters can be adjusted in the script:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `IDENTITY_THRESHOLD` | 0.90 | Min identity for conserved regions |
| `MIN_REGION_LEN` | 16 | Min conserved region length (bp) |
| `MAX_PRIMER_LEN` | 24 | Max primer length (bp) |
| `MAX_AMB` | 4 | Max ambiguous positions per primer |
| `MIN_TM` | 51.4 | Min melting temperature (°C) |
| `MIN_AMP_SEQ` | 500 | Min amplicon for sequencing (bp) |
| `MAX_AMP_SEQ` | 1200 | Max amplicon for sequencing (bp) |
| `MIN_AMP_QPCR` | 80 | Min amplicon for qPCR (bp) |
| `MAX_AMP_QPCR` | 350 | Max amplicon for qPCR (bp) |
| `SENS_TARGET` | 80.0 | Target true PCR sensitivity (%) |
| `COMBO_DEGEN_PREFER` | 24 | Preferred max combo degeneracy |
| `COMBO_DEGEN_MAX` | 32 | Hard max combo degeneracy |

## Output Format

### Primer Pool Table

| Location | Primer | Tm | VN% | Asia% | Global% |
|----------|--------|-----|-----|-------|---------|
| 65-207 | `TGCATTGGYTAYCATGCA` | 53.8 | 98.1 | 91.4 | 97.1 |

### Selected Primers

For each recommended primer pair:
- Forward primer sequence (use as-is for PCR)
- Reverse primer sequence (FASTA orientation)
- Reverse complement of reverse primer (order this for PCR)
- Amplicon size
- True PCR sensitivity (% sequences where BOTH primers match)

## Important Notes

### Primer Orientation

- All primers in output are in **5'→3' FASTA reference orientation**
- **Forward primers:** Use as-is for PCR ordering
- **Reverse primers:** Must be **reverse complemented** before ordering

### True PCR Sensitivity

Unlike individual primer sensitivity, true PCR sensitivity requires **both** forward and reverse primers to match the **same** sequence. This is the actual expected amplification rate.

### Degeneracy

Primer degeneracy = product of ambiguity at each position:
- R, Y, W, S, M, K = 2-fold
- B, D, H, V = 3-fold  
- N = 4-fold

Example: `ATGRYC` = 1×1×1×2×2×1 = 4-fold degeneracy

## Example Usage

### Influenza HA Gene

```python
GENE_NAME  = 'H5_VN'
ALN_FILE   = 'H5_VN_aligned.fasta'
CROSS_PANELS = {
    'Asia':   'H5_Asia.fasta',
    'Global': 'H5_global.fasta',
}
```

### Custom Thresholds

For highly variable genes, you may need to:
```python
IDENTITY_THRESHOLD = 0.80  # Lower threshold for more regions
MIN_PRIMER_SENS = 50.0     # Accept lower individual sensitivity
SENS_TARGET = 70.0         # Lower true PCR target
```

## Troubleshooting

### "No valid combo found"

- Lower `IDENTITY_THRESHOLD` to find more conserved regions
- Lower `MIN_PRIMER_SENS` to include more primers in selection
- Check if alignment quality is sufficient

### Low sensitivity on cross-panels

- Cross-panels may contain divergent strains
- Consider designing separate primers for different clades
- Use more representative sequences in design panel

### Primers failing Tm filter

- Regions may be AT-rich
- Try lowering `MIN_TM` slightly (not below 50°C)
- Increase `MAX_PRIMER_LEN` to allow longer primers

## Citation

If you use this pipeline in your research, please cite:

> [Your citation information here]

## License

[Your license information here]

## Contact

[Your contact information here]
