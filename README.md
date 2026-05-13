# AI4VPrimer

**Automated degenerate primer design for viral surveillance — with a no-code web interface.**

---

## What is AI4VPrimer?

AI4VPrimer is a computational pipeline that designs degenerate PCR primers from aligned viral sequences. It identifies conserved genomic regions, constructs primers with IUPAC ambiguity codes to capture sequence diversity, and selects optimal primer pairs for:

- **Full-length sequencing** (500–1200 bp overlapping amplicons)
- **qPCR detection** (80–350 bp amplicons)

The pipeline evaluates primers by melting temperature, homodimer stability, GC content, degeneracy, and cross-panel sensitivity — then recommends the best combinations automatically.

---

## Why Use the Web Interface?

Traditionally, updating primer design parameters requires opening Python source code, understanding variable names, and editing values manually. This creates barriers for lab scientists and collaborators who need results quickly.

**The AI4VPrimer web interface solves this by:**

- Letting users configure everything through a simple form — no coding required
- Clearly labeling each parameter with its meaning and units
- Marking required fields (`*`) vs. optional ones
- Saving settings directly into the pipeline with one click
- Running the pipeline and displaying results in the browser

**This means:**
- A new virus target? Just paste your aligned FASTA path and click Run.
- Need tighter Tm constraints? Change one number and re-run in seconds.
- Collaborator in another lab? They clone the repo, run one command, and get the interface.

No Python knowledge needed. No risk of accidentally breaking the code.

---

## How It Works

```
┌─────────────────────────────────────────────────────┐
│  Web Interface (browser)                            │
│  - Fill in gene name, file paths                    │
│  - Adjust Tm, amplicon size, sensitivity targets    │
│  - Click Save → Click Run                           │
└──────────────────────┬──────────────────────────────┘
                       │ writes parameters
                       ▼
┌─────────────────────────────────────────────────────┐
│  primer_pipeline.py                                 │
│  1. Load aligned sequences                          │
│  2. Find conserved regions (≥90% identity)          │
│  3. Build degenerate primers (IUPAC codes)          │
│  4. Filter: Tm, homodimer, GC, degeneracy           │
│  5. Select best pairs for sequencing & qPCR         │
│  6. Validate on cross-panels (optional)             │
│  7. Output markdown report                          │
└─────────────────────────────────────────────────────┘
```

---

## Quick Start

### 1. Install

Option 1- with git installed
```bash
git clone https://github.com/YOUR_USERNAME/AI4VPrimer.git
cd AI4VPrimer
pip install -r requirements.txt
```
Option 2 — Download ZIP (easier for non-programmers):

   1. Go to your GitHub repo page
   2. Click the green Code button
   3. Click Download ZIP
   4. Unzip the folder

### 2. Launch the Interface

```bash
python app.py
```

The browser opens automatically to the AI4VPrimer configuration page.

> If it doesn't open, go to **http://127.0.0.1:5000** manually.

### 3. Configure

Fill in the required fields:

| Field | Description |
|-------|-------------|
| GENE_NAME * | Target name (e.g., `H5_VN`, `N1_Asia`) |
| ALN_FILE * | Path to your aligned FASTA file |
| OUT_FILE * | Path for the output report |
| PRIMERS_PER_REGION * | Primers to keep per conserved region (default: 3) |
| CROSS_PANELS | Optional: sensitivity check panels (JSON format) |
| SPECIFICITY_FILE | Optional: non-target sequences for specificity |

Adjust design parameters (Tm, amplicon sizes, etc.) or leave defaults.

### 4. Save & Run

1. Click **💾 Save Parameters**
2. Click **▶️ Run Pipeline**
3. View results in the browser or open the output file

### 5. Stop

Press **Ctrl+C** in Terminal when finished.

---

## Requirements

- Python 3.7+
- Flask ≥ 3.0
- Biopython
- NumPy
- primer3-py

Install all with:

```bash
pip install -r requirements.txt
```

---

## Input Requirements

- **Aligned FASTA file** (required): Multiple sequences aligned with MAFFT, MUSCLE, or similar. The pipeline uses the alignment to identify conserved regions.
- **Cross-panel FASTA files** (optional): Unaligned sequences from broader geographic/temporal panels for sensitivity validation.
- **Specificity FASTA file** (optional): Non-target sequences to verify primers don't cross-react.

---

## Output

The pipeline generates a markdown report containing:

- **Primer Pool**: All primers passing filters with location, Tm, GC%, sensitivity
- **Sequencing Selection**: Two overlapping primer pairs covering the target region
- **qPCR Selection**: Single best primer pair for detection assays
- **Reverse complements**: Ready for PCR ordering

---

## File Structure

```
AI4VPrimer/
├── app.py                 # Web server (launches the interface)
├── index.html             # Web interface page
├── primer_pipeline.py     # Core pipeline code
├── requirements.txt       # Python dependencies
└── README.md              # This file
```

---

## License

MIT
