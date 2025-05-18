# 📘 Primer Analysis and FASTA Management Notebooks

This repository contains a set of Jupyter notebooks designed to assist with primer design, analysis, and FASTA file management in bioinformatics workflows.

## 📂 Notebook Summaries

---

### 1. `In_silicon_aijia_degeneracy.ipynb`

**Purpose**:  
Expands degenerate primers into all possible non-degenerate sequences and evaluates their effectiveness.

**Key Features**:
- `expand_degenerate_primer(primer)`: Generates all non-degenerate combinations from a degenerate primer.
- `expand(base)`: Maps individual degenerate bases to their non-degenerate equivalents.
- `do_primers_work(FP)`: Evaluates a forward primer’s effectiveness based on sequence analysis.

**Use Case**:  
Ideal for researchers working with degenerate primers to ensure proper expansion and testing before biological application.

---

### 2. `In_silicon_aijia.ipynb`

**Purpose**:  
Checks the **sensitivity** of degenerate primer pairs.

---

### 3. `sensitivity and specificity.ipynb`

**Purpose**:  
Calculates sensitivity and specificity of primer matches against sequence data.

**Key Features**:
- `calculate_metrics(fasta_file1, FP, RP, merged_fasta_file)`: Computes sensitivity and specificity metrics using provided primers and merged FASTA files.

**Use Case**:  
Useful for assessing the **precision and recall** of primer-based detection—important for making data-driven decisions about primer performance.

---

### 4. `merge fasta.ipynb`

**Purpose**:  
Merges multiple FASTA files into a single consolidated file for easier data handling.

**Key Features**:
- `merge_fasta_files(file_list, output_file)`: Combines FASTA sequences from several files into one output file.

**Use Case**:  
Great for managing large datasets by simplifying FASTA file consolidation for further analysis or sharing.

---
