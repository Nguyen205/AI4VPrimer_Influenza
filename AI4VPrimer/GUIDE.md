# AI4VPrimer — Detailed Usage Guide

---

## 1. Installation

### Option A: Git clone
```bash
git clone https://github.com/YOUR_USERNAME/AI4VPrimer.git
cd AI4VPrimer
pip install -r requirements.txt
```

### Option B: Download ZIP
1. Go to the GitHub repository page
2. Click green **Code** button → **Download ZIP**
3. Unzip the folder
4. Open Terminal, `cd` into the folder
5. Run `pip install -r requirements.txt`

> ⚠️ Flask 3.0 or higher is required. Check: `python -c "import flask; print(flask.__version__)"`

---

## 2. Launching the Web Interface

```bash
cd AI4VPrimer
python app.py
```

The browser opens automatically to `http://127.0.0.1:5000`.

> If it doesn't open, manually navigate to **http://127.0.0.1:5000**

To stop the server when finished: press **Ctrl+C** in Terminal.

---

## 3. Web Interface Walkthrough

The interface has two sections: **Configuration** and **Design & Selection Parameters**.

### Configuration (required fields marked with *)

| Field | What to enter | Example |
|-------|---------------|---------|
| GENE_NAME * | Name for your target | `H5_VN` |
| ALN_FILE * | Full path to aligned FASTA | `/Users/you/data/H5_aligned.fasta` |
| OUT_FILE * | Full path for output (must end in `.md`) | `/Users/you/Desktop/H5_output.md` |
| PRIMERS_PER_REGION * | How many primers to keep per conserved region | `3` |
| SUBSAMPLE_REPS | Repetitions for large datasets (0 = auto) | `0` |
| CROSS_PANELS | Optional: extra panels for sensitivity check | `{"Asia": "/path/to/asia.fasta"}` |
| SPECIFICITY_FILE | Optional: non-target sequences | `/path/to/non_targets.fasta` |

**Important:**
- All file paths must be **full/absolute paths** (e.g., `/Users/you/...`)
- ALN_FILE must be **pre-aligned** (all sequences same length)
- OUT_FILE must end with **`.md`**
- CROSS_PANELS uses JSON format: `{"name": "/full/path.fasta", "name2": "/path2.fasta"}`
- Leave CROSS_PANELS and SPECIFICITY_FILE empty if not needed

### Design & Selection Parameters

| Parameter | Default | What it controls |
|-----------|---------|------------------|
| MAX_PRIMER_LEN | 24 | Maximum primer length in bases |
| MIN_TM | 51.4 | Minimum melting temperature (°C). Primers below this are rejected |
| MAX_LOW_IDENTITY_BRIDGE | 1 | How many low-identity positions can be bridged within a conserved region |
| MIN_AMP_SEQ | 500 | Minimum amplicon size for sequencing pairs (bp) |
| MAX_AMP_SEQ | 1200 | Maximum amplicon size for sequencing pairs (bp) |
| MIN_AMP_QPCR | 80 | Minimum amplicon size for qPCR pair (bp) |
| MAX_AMP_QPCR | 350 | Maximum amplicon size for qPCR pair (bp) |
| MIN_PRIMER_SENS | 60 | Minimum individual primer sensitivity (%) to enter the pool |
| SENS_TARGET | 80 | Target true PCR sensitivity (%) for pair selection |
| OVERLAP_LIMITS | [200,300,400,500] | Progressive overlap relaxation for sequencing pairs (bp) |
| COMBO_DEGEN_PREFER | 24 | Preferred max combined degeneracy (forward × reverse) |

---

## 4. Running the Pipeline

1. Fill in all required fields
2. Click **💾 Save Parameters** — writes your settings into the pipeline code
3. Click **▶️ Run Pipeline** — executes the pipeline

> ⚠️ You must Save before Run. Otherwise the pipeline uses old/empty values.

### What happens during a run:

1. **Load sequences** — reads your aligned FASTA
2. **Check dataset size:**
   - ≤1000 sequences → standard mode (designs on all sequences)
   - &gt;1000 sequences → subsample mode (designs on 1000-sequence subsets, repeats N times)
3. **Find conserved regions** — scans alignment for stretches with ≥90% identity
4. **Design primers** — builds degenerate primers, trims to optimal length, checks Tm and dimer
5. **Validate** — (subsample mode) tests all primers against the full dataset
6. **Select pairs** — picks best combinations for sequencing and qPCR
7. **Write report** — saves output to your specified path

---

## 5. Understanding the Output

The output `.md` file contains:

### Primer Pool Table

| Column | Meaning |
|--------|---------|
| Location | Position in the alignment (start-end) |
| Primer | The primer sequence (with IUPAC codes) |
| Tm | Melting temperature (°C) |
| GC% | GC content percentage |
| VN% | Sensitivity — % of sequences matched |

### Sequencing Selection

Two overlapping primer pairs that together cover the maximum span of the gene:
- **C1 Fwd/Rev** — first amplicon (500-1200 bp)
- **C2 Fwd/Rev** — second amplicon (500-1200 bp, overlapping with C1)
- **True PCR %** — % of sequences where BOTH primers in the pair match

### qPCR Selection

One primer pair producing a short amplicon (80-350 bp) with the highest sensitivity.

### Reverse Complements

The output includes RC sequences for reverse primers — these are what you order for PCR.

---

## 6. Parameter Tuning Guide

### When to adjust MIN_TM:
- **Increase** (e.g., 55°C) if you need more stringent annealing, fewer off-target products
- **Decrease** (e.g., 50°C) if the pipeline finds too few primers (highly variable target)

### When to adjust amplicon sizes:
- **Sequencing (MIN/MAX_AMP_SEQ):** Match your sequencing platform's read length
- **qPCR (MIN/MAX_AMP_QPCR):** Standard qPCR works best at 80-200 bp

### When to adjust PRIMERS_PER_REGION:
- **1** — fastest, picks only the best primer per region
- **3** — gives more options for pair selection (recommended)
- **5+** — for very diverse targets where you want maximum coverage

### When to adjust SUBSAMPLE_REPS:
- **0 (auto)** — uses ceil(total/1000) repetitions
- **Manual (2, 3, 5)** — more reps = more primer candidates but longer runtime
- Each rep takes ~5 minutes for a 1500 bp gene

### When to adjust COMBO_DEGEN_PREFER:
- **Lower (8-16)** — fewer primer variants in the tube, cleaner PCR but may miss some sequences
- **Higher (24-32)** — more variants, higher sensitivity but potentially more non-specific products

---

## 7. Subsample Mode (>1000 sequences)

When your dataset exceeds 1000 sequences, the pipeline automatically:

1. Randomly picks 1000 sequences → designs primers
2. Repeats with different random subsets (N times)
3. Pools all unique primers from all repetitions
4. Tests every pooled primer against ALL sequences for true sensitivity
5. Keeps only primers with sensitivity ≥ MIN_PRIMER_SENS
6. Selects best pairs from the validated pool

This keeps runtime manageable (~5 min per repetition) while ensuring primers are validated on the full dataset.

---

## 8. Troubleshooting

| Error | Cause | Fix |
|-------|-------|-----|
| ALN_FILE not found | Wrong path | Use full absolute path |
| Sequences have different lengths | File not aligned | Align with MAFFT/MUSCLE first |
| OUT_FILE must end with .md | Wrong extension | Add `.md` to output path |
| No valid qPCR pair found | No primers fit amplicon range | Try increasing MAX_AMP_QPCR or decreasing MIN_TM |
| No valid overlapping combo | Primers too far apart or too close | Try adjusting MIN/MAX_AMP_SEQ |
| Pipeline timed out | Too many sequences in standard mode | Should auto-switch to subsample mode; check SUBSAMPLE_REPS |
| slice indices error | PRIMERS_PER_REGION saved as text | Re-enter as a number and Save again |

---

## 9. IUPAC Ambiguity Codes Reference

| Code | Bases | Meaning |
|------|-------|---------|
| R | A, G | Purine |
| Y | C, T | Pyrimidine |
| W | A, T | Weak |
| S | C, G | Strong |
| M | A, C | Amino |
| K | G, T | Keto |
| B | C, G, T | Not A |
| D | A, G, T | Not C |
| H | A, C, T | Not G |
| V | A, C, G | Not T |
| N | A, C, G, T | Any |

Degeneracy = number of actual sequences a degenerate primer represents.
Example: `ATCRG` (R=2 variants) → degeneracy = 2 → represents `ATCAG` and `ATCGG`.
