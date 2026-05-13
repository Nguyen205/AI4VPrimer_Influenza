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
pip install -r requirements.txt
```

Or install manually:

```bash
pip install biopython numpy primer3-py "flask>=3.0"
```

> ⚠️ **Flask 3.0 or higher is required.** Older versions (e.g., 1.x, 2.x) will fail due to dependency incompatibilities. Check your version with `python -c "import flask; print(flask.__version__)"`

## Files

| File | Purpose |
|------|---------|
| `primer_pipeline.py` | Main design & selection pipeline |
| `app.py` | Web API for editing pipeline parameters |
| `index.html` | Browser-based parameter editor |
| `sensitivity and specificity.ipynb` | Additional validation of primers |
| `merge fasta.ipynb` | Combine multiple FASTA files |

---

## Web Interface (Parameter Editor)

A browser-based interface lets anyone edit pipeline parameters without reading the code.

### Quick Start

```bash
pip install flask
python app.py
```

Then open **http://localhost:5000** in your browser.

### How It Works

- The web page shows all CONFIG and PARAMETERS fields with descriptions
- Edit values in the form and click **Save**
- The API updates `primer_pipeline.py` directly (the actual source code)
- Then run `python primer_pipeline.py` as usual with the new settings

### API Endpoints

| Method | Endpoint | Description |
|--------|----------|-------------|
| GET | `/api/params` | Returns current CONFIG and PARAMETERS as JSON |
| POST | `/api/params` | Updates variables in primer_pipeline.py |

**POST example:**
```bash
curl -X POST http://localhost:5000/api/params \
  -H "Content-Type: application/json" \
  -d '{"MIN_TM": 55.0, "GENE_NAME": "N1_VN"}'
```

---

## Getting Started (Step-by-Step)

### Step 1: Clone the Repository

```bash
git clone https://github.com/YOUR_USERNAME/pipeline-for-primer-design.git
cd pipeline-for-primer-design
```

### Step 2: Install Dependencies

```bash
pip install -r requirements.txt
```

> ⚠️ **Flask 3.0 or higher is required.** Check with: `python -c "import flask; print(flask.__version__)"`

### Step 3: Prepare Your Input File

- Align your target sequences using MAFFT, MUSCLE, or another aligner
- Save the aligned output as a `.fasta` file on your computer
- Note the full path to this file (e.g., `/Users/you/data/my_aligned.fasta`)

### Step 4: Launch the Web Interface

Open **Terminal** and type:

```bash
cd pipeline-for-primer-design
python app.py
```

A browser window will open automatically at `http://127.0.0.1:5000` showing the **AI4VPrimer** page.

> If the browser doesn't open automatically, manually go to **http://127.0.0.1:5000**

### Step 5: Fill in the Configuration

On the webpage, fill in the required fields (marked with `*`):

| Field | What to enter |
|-------|---------------|
| GENE_NAME * | Your target name (e.g., `H5_VN`) |
| ALN_FILE * | Full path to your aligned FASTA file |
| OUT_FILE * | Full path where you want the output report saved |
| PRIMERS_PER_REGION * | Number of primers per region (default: 3) |
| CROSS_PANELS | Optional: JSON format, e.g. `{"Asia": "/path/to/asia.fasta"}` |
| SPECIFICITY_FILE | Optional: path to non-target FASTA for specificity check |

Adjust the **Design & Selection Parameters** below if needed, or leave the defaults.

### Step 6: Save and Run

1. Click **💾 Save Parameters** — this writes your settings into the pipeline code
2. Click **▶️ Run Pipeline** — this executes the pipeline and shows results in the browser

> ⚠️ You must click **Save** before **Run**. Otherwise the pipeline runs with empty/old values.

### Step 7: View Results

The output report is saved to the path you specified in `OUT_FILE`. It contains:
- A table of all candidate primers with Tm, GC%, and sensitivity
- Recommended primer pairs for full-length sequencing
- Recommended primer pair for qPCR
- Reverse complements ready for PCR ordering

### Step 8: Stop the Server

When finished, go back to Terminal and press **Ctrl+C** to stop the web server.

---

## Advanced Usage

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
