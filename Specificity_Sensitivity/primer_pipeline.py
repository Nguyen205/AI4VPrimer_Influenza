"""
Unified Primer Design & Selection Pipeline
==========================================
A single pipeline that:
  1. Designs degenerate primers from alignment (relaxation approach)
  2. Filters by Tm and homodimer
  3. Selects optimal primer pairs for:
     - Full-length sequencing (500-1200 bp amplicons, overlapping)
     - qPCR (80-350 bp amplicons)

Usage:
  Modify the CONFIG section below for your gene/files, then run.
"""

import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
import re, itertools

# ══════════════════════════════════════════════════════════════════════════════
# CONFIG - Modify this section for your gene
# ══════════════════════════════════════════════════════════════════════════════
GENE_NAME  = 'H5_VN'
ALN_FILE   = '/Users/aijiazhou/Desktop/Markdown/H5_VN_aligned.fasta'

# Cross-panels for sensitivity check (set to None or {} if not needed)
CROSS_PANELS = {
    'Asia':   '/Users/aijiazhou/Desktop/Fasta_file/H5_Asia (1).fasta',
    'Global': '/Users/aijiazhou/Desktop/Fasta_file/H5_global.fasta',
}

# Specificity check file (set to None if not needed)
# Specificity = % of non-target sequences NOT matched by primer combo
SPECIFICITY_FILE = None  # e.g., '/path/to/non_target_sequences.fasta'

OUT_FILE   = '/Users/aijiazhou/Desktop/Pipeline for primer design/primer_output.md'

# Number of primers to keep per conserved region (if enough pass filters)
PRIMERS_PER_REGION = 3

# ══════════════════════════════════════════════════════════════════════════════
# PARAMETERS
# ══════════════════════════════════════════════════════════════════════════════
# Design parameters
IDENTITY_THRESHOLD      = 0.90
MIN_REGION_LEN          = 16
MAX_PRIMER_LEN          = 24
MAX_AMB                 = 4
IMPROVE_THRESH          = 0.02
MIN_TM                  = 51.4
MAX_LOW_IDENTITY_BRIDGE = 1

# Selection parameters
MIN_AMP_SEQ      = 500    # min amplicon for sequencing
MAX_AMP_SEQ      = 1200   # max amplicon for sequencing
MIN_AMP_QPCR     = 80     # min amplicon for qPCR
MAX_AMP_QPCR     = 350    # max amplicon for qPCR
MIN_PRIMER_SENS  = 60.0   # min individual primer sensitivity to enter pool
SENS_TARGET      = 80.0   # target true PCR sensitivity
OVERLAP_LIMITS   = [200, 300, 400, 500]  # progressive relaxation
COMBO_DEGEN_PREFER = 24   # preferred max combo degeneracy (fwd × rev)
COMBO_DEGEN_MAX    = 32   # hard max combo degeneracy

# ══════════════════════════════════════════════════════════════════════════════
# IUPAC & HELPERS
# ══════════════════════════════════════════════════════════════════════════════
IUPAC = {'R':'[AG]','Y':'[CT]','W':'[AT]','S':'[CG]','M':'[AC]','K':'[GT]',
         'B':'[CGT]','D':'[AGT]','H':'[ACT]','V':'[ACG]','N':'[ACGT]'}
IUPAC_MAP = {frozenset(['A','G']):'R', frozenset(['C','T']):'Y',
             frozenset(['A','T']):'W', frozenset(['C','G']):'S',
             frozenset(['A','C']):'M', frozenset(['G','T']):'K',
             frozenset(['C','G','T']):'B', frozenset(['A','G','T']):'D',
             frozenset(['A','C','T']):'H', frozenset(['A','C','G']):'V',
             frozenset(['A','C','G','T']):'N'}
AMBIG_CODES = ['R','Y','W','S','M','K','B','D','H','V']
UNAMBIG = set('ATCG')
DEGEN_VAL = {'R':2,'Y':2,'W':2,'S':2,'M':2,'K':2,'B':3,'D':3,'H':3,'V':3,'N':4,'A':1,'T':1,'C':1,'G':1}
nucs = ['A','T','C','G']

def p2re(p): return ''.join(IUPAC.get(c, c) for c in p.upper())
def rc(s): return str(Seq(s).reverse_complement())

class _Rec:
    def __init__(self, seq): self.seq = seq

def calc_hit(primer, seqs):
    fwd = re.compile(p2re(primer), re.I)
    rev = re.compile(p2re(rc(primer)), re.I)
    return sum(1 for s in seqs if fwd.search(str(s.seq)) or rev.search(str(s.seq)))

def calc_tm(seq):
    IUPAC_MIN = {'R':'A','Y':'T','W':'A','S':'C','M':'A','K':'T','B':'T','D':'A','H':'A','V':'A','N':'A'}
    IUPAC_MAX = {'R':'G','Y':'C','W':'T','S':'G','M':'C','K':'G','B':'G','D':'G','H':'C','V':'G','N':'G'}
    s_min = ''.join(IUPAC_MIN.get(c.upper(), c.upper()) for c in seq)
    s_max = ''.join(IUPAC_MAX.get(c.upper(), c.upper()) for c in seq)
    try:
        return round((mt.Tm_NN(s_min, Na=50, dnac1=300, dnac2=0) + mt.Tm_NN(s_max, Na=50, dnac1=300, dnac2=0)) / 2, 1)
    except:
        return round(mt.Tm_Wallace(s_min), 1)

def dimer_ok(primer):
    from primer3 import calc_homodimer
    def deambig(s):
        m={'R':'A','Y':'C','W':'A','S':'C','M':'A','K':'G','B':'C','D':'A','H':'A','V':'A','N':'A'}
        return ''.join(m.get(c.upper(),c.upper()) for c in s)
    for seq in [primer, rc(primer)]:
        r = calc_homodimer(deambig(seq), mv_conc=100, dv_conc=2, dntp_conc=0.8, dna_conc=15)
        if r.tm > 35 or r.dg < -9000: return False
    return True

def ends_clean(p, n=2):
    return sum(1 for c in p[:n] if c.upper() not in UNAMBIG) == 0 and sum(1 for c in p[-n:] if c.upper() not in UNAMBIG) == 0

def primer_passes(primer):
    return calc_tm(primer) >= MIN_TM and calc_tm(rc(primer)) >= MIN_TM and dimer_ok(primer)

def calc_gc_content(primer):
    """Calculate GC content (0-1). For ambiguous bases, use average."""
    gc_map = {'G':1, 'C':1, 'A':0, 'T':0, 'S':1, 'R':0.5, 'Y':0.5, 'W':0, 'M':0.5, 'K':0.5,
              'B':0.67, 'D':0.33, 'H':0.33, 'V':0.67, 'N':0.5}
    return sum(gc_map.get(c.upper(), 0.5) for c in primer) / len(primer)

def has_consecutive_bases(primer, n=4):
    """Check if primer has n or more consecutive identical bases (AAAA, CCCC, etc.)."""
    for base in 'ATCG':
        if base * n in primer.upper():
            return True
    return False

def primer_quality_score(primer):
    """Score primer quality for tiebreaking. Higher is better."""
    gc = calc_gc_content(primer)
    # Ideal GC is 40-60% (0.4-0.6), penalize deviation
    gc_score = 1 - abs(gc - 0.5) * 2  # 1.0 at 50%, 0.8 at 40% or 60%, etc.
    # Penalize consecutive bases
    consec_penalty = -0.5 if has_consecutive_bases(primer, 4) else 0
    return gc_score + consec_penalty

def true_pcr_sens(fseq, rseq, seqs):
    ff = re.compile(p2re(fseq), re.I); fr = re.compile(p2re(rc(fseq)), re.I)
    rf = re.compile(p2re(rseq), re.I); rr = re.compile(p2re(rc(rseq)), re.I)
    hits = sum(1 for s in seqs if (ff.search(s) or fr.search(s)) and (rf.search(s) or rr.search(s)))
    return round(100 * hits / len(seqs), 1)

# ══════════════════════════════════════════════════════════════════════════════
# PRIMER DESIGN FUNCTIONS
# ══════════════════════════════════════════════════════════════════════════════
def build_primer(aln, start, end):
    """Build primer with IUPAC codes for positions where minor base >=5%."""
    primer, amb = [], 0
    for i in range(start, end):
        col = [s[i] for s in aln if s[i] != '-']
        if not col: continue
        counts = {n: col.count(n) for n in nucs}
        t = sum(counts.values())
        dominant = max(counts, key=counts.get)
        if counts[dominant] / t >= IDENTITY_THRESHOLD:
            primer.append(dominant)
        else:
            # Include bases with >=5% frequency
            present = frozenset(n for n in nucs if counts[n]/t >= 0.05)
            if len(present) <= 1: 
                primer.append(dominant)
            elif amb < MAX_AMB: 
                primer.append(IUPAC_MAP.get(present, 'N'))
                amb += 1
            else: 
                primer.append(dominant)
    return ''.join(primer)

def best_trim(primer, seqs, total):
    best_p, best_hit, fallback_p, fallback_hit = None, -1, None, -1
    for length in range(MIN_REGION_LEN, MAX_PRIMER_LEN + 1):
        for start in range(len(primer) - length + 1):
            candidate = primer[start:start + length]
            hit = calc_hit(candidate, seqs)
            if calc_tm(candidate) >= MIN_TM and calc_tm(rc(candidate)) >= MIN_TM and hit > best_hit:
                best_hit, best_p = hit, candidate
            if hit > fallback_hit: fallback_hit, fallback_p = hit, candidate
    return (best_p, best_hit) if best_p else (fallback_p or primer[:MAX_PRIMER_LEN], fallback_hit if fallback_p else calc_hit(primer[:MAX_PRIMER_LEN], seqs))

def reduce_degeneracy(primer, base_hit, seqs, total, aln, aln_start, aln_len):
    """Remove ambiguous bases if sensitivity loss <2%."""
    result, cur_hit, changed = list(primer), base_hit, True
    while changed:
        changed = False
        for i in range(len(result)):
            if result[i].upper() in UNAMBIG: continue
            col = [s[aln_start + i] for s in aln if aln_start + i < len(s) and s[aln_start + i] != '-']
            dom = max('ATCG', key=lambda n: col.count(n)) if col else result[i]
            cand = result[:]; cand[i] = dom
            hit = calc_hit(''.join(cand), seqs)
            if (cur_hit - hit) / total < 0.02: result[i] = dom; cur_hit = hit; changed = True; break
    def perms(p): return np.prod([DEGEN_VAL.get(c.upper(), 1) for c in p])
    while perms(result) > 24:
        amb = [(i, c) for i, c in enumerate(result) if c.upper() not in UNAMBIG]
        if not amb: break
        i, c = max(amb, key=lambda x: DEGEN_VAL.get(x[1].upper(), 1))
        col = [s[aln_start + i] for s in aln if aln_start + i < len(s) and s[aln_start + i] != '-']
        result[i] = max('ATCG', key=lambda n: col.count(n)) if col else 'A'
    return ''.join(result), calc_hit(''.join(result), seqs)

def dominant_bases(primer, raw_seqs):
    fwd, rev = re.compile(p2re(primer), re.I), re.compile(p2re(rc(primer)), re.I)
    cols = [[] for _ in range(len(primer))]
    for s in raw_seqs:
        seq = str(s.seq).upper()
        m = fwd.search(seq)
        if m:
            for i, b in enumerate(m.group()):
                if b in UNAMBIG: cols[i].append(b)
        else:
            m = rev.search(seq)
            if m:
                matched = str(Seq(m.group()).reverse_complement()).upper()
                for i, b in enumerate(matched):
                    if b in UNAMBIG: cols[i].append(b)
    return [max(set(c), key=c.count) if c else 'A' for c in cols]

def fix_ends_replace(primer, raw_seqs, n=2):
    dom, p = dominant_bases(primer, raw_seqs), list(primer)
    for i in list(range(n)) + list(range(len(p)-n, len(p))):
        if p[i].upper() not in UNAMBIG: p[i] = dom[i]
    return ''.join(p)

def get_consensus_base(aln, pos):
    col = [s[pos] for s in aln if s[pos] != '-']
    return max(nucs, key=lambda n: col.count(n)) if col else 'A'

def find_aln_pos(primer, aln, reg_start, reg_end, aln_len):
    window = max(0, reg_start - 50)
    cons = ''.join(get_consensus_base(aln, i) for i in range(window, min(aln_len, reg_end + 50)))
    m = re.search(re.sub('[^ATCG]', '.', primer.upper()), cons)
    return (window + m.start(), window + m.end()) if m else (reg_start, reg_start + len(primer))

def fix_ends_add(primer, aln_start, aln_end, aln, aln_len, n=2):
    p, s, e = primer, aln_start, aln_end
    for _ in range(10):
        if ends_clean(p, n): break
        if any(c.upper() not in UNAMBIG for c in p[-n:]) and e < aln_len: p = p + get_consensus_base(aln, e); e += 1
        if any(c.upper() not in UNAMBIG for c in p[:n]) and s > 0: s -= 1; p = get_consensus_base(aln, s) + p
    return p

# ══════════════════════════════════════════════════════════════════════════════
# PRIMER SELECTION FUNCTIONS
# ══════════════════════════════════════════════════════════════════════════════
def calc_degeneracy(primer):
    """Calculate total degeneracy of a primer."""
    return np.prod([DEGEN_VAL.get(c.upper(), 1) for c in primer])

def combo_degeneracy(fwd, rev):
    """Calculate combined degeneracy of a primer pair."""
    return calc_degeneracy(fwd) * calc_degeneracy(rev)

def build_combos(primers, min_amp, max_amp):
    combos = []
    for f in primers:
        for r in primers:
            if f is r or f['pos_start'] >= r['pos_start']: continue
            amp = r['pos_end'] - f['pos_start']
            if not (min_amp <= amp <= max_amp): continue
            if f['vn_pct'] < MIN_PRIMER_SENS or r['vn_pct'] < MIN_PRIMER_SENS: continue
            # Filter by combo degeneracy
            degen = combo_degeneracy(f['primer'], r['primer'])
            if degen > COMBO_DEGEN_MAX: continue
            combos.append({'f': f, 'r': r, 'amp': amp, 'degen': degen})
    return combos

def find_overlap_pairs(combos, max_overlap):
    pairs = []
    for a, b in itertools.combinations(combos, 2):
        if a['f']['pos_start'] > b['f']['pos_start']: a, b = b, a
        if not (a['f']['pos_start'] < b['f']['pos_start'] < a['r']['pos_end'] < b['r']['pos_end']): continue
        if a['f']['primer'] == b['f']['primer'] or a['r']['primer'] == b['r']['primer']: continue
        if a['f']['primer'] == b['r']['primer'] or a['r']['primer'] == b['f']['primer']: continue
        overlap = a['r']['pos_end'] - b['f']['pos_start']
        if overlap > max_overlap: continue
        pairs.append((a, b, overlap, b['r']['pos_end'] - a['f']['pos_start']))
    return sorted(pairs, key=lambda x: -x[3])

def select_sequencing(primers, vn_seqs):
    best_result = None
    for limit in OVERLAP_LIMITS:
        combos = build_combos(primers, MIN_AMP_SEQ, MAX_AMP_SEQ)
        pairs = find_overlap_pairs(combos, limit)
        for a, b, overlap, span in pairs[:20]:
            # Prefer combos with degeneracy ≤24, but allow up to 32
            total_degen = a['degen'] + b['degen']
            s1 = true_pcr_sens(a['f']['primer'], a['r']['primer'], vn_seqs)
            s2 = true_pcr_sens(b['f']['primer'], b['r']['primer'], vn_seqs)
            # Score: prioritize sensitivity, then prefer lower degeneracy
            score = min(s1, s2) - (0.1 if total_degen > COMBO_DEGEN_PREFER * 2 else 0)
            if best_result is None or score > min(best_result[4], best_result[5]) - (0.1 if (best_result[0]['degen'] + best_result[1]['degen']) > COMBO_DEGEN_PREFER * 2 else 0):
                best_result = (a, b, overlap, span, s1, s2, limit)
            if s1 >= SENS_TARGET and s2 >= SENS_TARGET and a['degen'] <= COMBO_DEGEN_PREFER and b['degen'] <= COMBO_DEGEN_PREFER:
                return best_result
    return best_result

def select_qpcr(primers, vn_seqs):
    combos = build_combos(primers, MIN_AMP_QPCR, MAX_AMP_QPCR)
    best, best_score = None, -1
    for c in combos:
        sens = true_pcr_sens(c['f']['primer'], c['r']['primer'], vn_seqs)
        # Prefer degeneracy ≤24, penalize higher
        score = sens - (1 if c['degen'] > COMBO_DEGEN_PREFER else 0)
        if score > best_score: best, best_score = (c, sens), score
    return best if best else None

# ══════════════════════════════════════════════════════════════════════════════
# MAIN PIPELINE
# ══════════════════════════════════════════════════════════════════════════════
if __name__ == '__main__':
    # ── Load data ─────────────────────────────────────────────────────────────
    print(f"=== {GENE_NAME} Primer Design & Selection Pipeline ===\n")
    print("Loading alignment...")
    _records = list(SeqIO.parse(ALN_FILE, 'fasta'))
    aln = [str(s.seq).upper() for s in _records]
    aln_len = len(aln[0])
    raw_seqs = [_Rec(s.replace('-', '')) for s in aln]
    total = len(raw_seqs)
    vn_seqs = [str(s.seq) for s in raw_seqs]
    print(f"  Design panel: {total} sequences, alignment: {aln_len} bp")

    # Load cross-panels (optional)
    cross_seqs = {}
    if CROSS_PANELS:
        for name, path in CROSS_PANELS.items():
            cross_seqs[name] = [str(r.seq).replace('-','').upper() for r in SeqIO.parse(path, 'fasta')]
            print(f"  {name} panel: {len(cross_seqs[name])} sequences")

    # Load specificity panel (optional)
    spec_seqs = None
    if SPECIFICITY_FILE:
        spec_seqs = [str(r.seq).replace('-','').upper() for r in SeqIO.parse(SPECIFICITY_FILE, 'fasta')]
        print(f"  Specificity panel: {len(spec_seqs)} sequences")

    # ── Compute identity ──────────────────────────────────────────────────────
    freq = {n: np.zeros(aln_len) for n in nucs + ['-']}
    for seq in aln:
        for i, base in enumerate(seq):
            if base in freq: freq[base][i] += 1
    for k in freq: freq[k] /= len(aln)
    identity = np.array([max(freq[n][i] for n in nucs) for i in range(aln_len)])

    # ── Extract regions ───────────────────────────────────────────────────────
    print(f"\nExtracting conserved regions (>={IDENTITY_THRESHOLD*100:.0f}% identity, >={MIN_REGION_LEN} bp)...")
    regions, i = [], 0
    while i < aln_len:
        if identity[i] >= IDENTITY_THRESHOLD:
            start, low_count, j = i, 0, i + 1
            while j < aln_len:
                if identity[j] >= IDENTITY_THRESHOLD: low_count = 0; j += 1
                elif low_count < MAX_LOW_IDENTITY_BRIDGE: low_count += 1; j += 1
                else: break
            end = j
            while end > start and identity[end - 1] < IDENTITY_THRESHOLD: end -= 1
            if (end - start) >= MIN_REGION_LEN: regions.append((start, end))
            i = j
        else: i += 1
    print(f"  Found {len(regions)} regions")

    # ── Design primers ────────────────────────────────────────────────────────
    print("\nDesigning primers...")
    primer_pool = []
    for reg_start, reg_end in regions:
        raw_primer = build_primer(aln, reg_start, reg_end)
        trimmed, _ = best_trim(raw_primer, raw_seqs, total)
        
        # Reduce degeneracy (remove ambiguous bases if sensitivity loss <2%)
        aln_s0, _ = find_aln_pos(trimmed, aln, reg_start, reg_end, aln_len)
        trimmed, _ = reduce_degeneracy(trimmed, calc_hit(trimmed, raw_seqs), raw_seqs, total, aln, aln_s0, aln_len)

        # Fix ends with two strategies
        rep = fix_ends_replace(trimmed, raw_seqs)
        aln_s, aln_e = find_aln_pos(trimmed, aln, reg_start, reg_end, aln_len)
        add = fix_ends_add(trimmed, aln_s, aln_e, aln, aln_len)

        # Collect candidates for this region
        region_candidates = []
        for p, strat in [(rep, 'replace'), (add, 'add-nt')]:
            if not primer_passes(p): continue
            hit = calc_hit(p, raw_seqs)
            tm = calc_tm(p)
            cross_pct = {name: 100 * calc_hit(p, [_Rec(s) for s in seqs]) / len(seqs) for name, seqs in cross_seqs.items()} if cross_seqs else {}
            gc = calc_gc_content(p)
            quality = primer_quality_score(p)
            region_candidates.append({
                'primer': p, 'tm': tm, 'pos_start': reg_start, 'pos_end': reg_end,
                'vn_pct': 100 * hit / total, 'cross_pct': cross_pct, 'strategy': strat,
                'gc': gc, 'quality': quality
            })
        
        # Sort by sensitivity (desc), then quality score (desc) for tiebreaking
        region_candidates.sort(key=lambda x: (-x['vn_pct'], -x['quality']))
        
        # Keep up to PRIMERS_PER_REGION primers from this region
        for p in region_candidates[:PRIMERS_PER_REGION]:
            primer_pool.append(p)
            print(f"  {reg_start}-{reg_end} ({p['strategy']}): {p['primer']} Tm={p['tm']} VN={p['vn_pct']:.1f}% GC={p['gc']*100:.0f}%")

    # Remove duplicates (keep first occurrence which has better score)
    seen = set()
    primer_pool = [p for p in primer_pool if not (p['primer'] in seen or seen.add(p['primer']))]
    print(f"\nPrimer pool: {len(primer_pool)} unique primers passing filters")

    # ── Selection ─────────────────────────────────────────────────────────────
    print("\n" + "="*80)
    print("PRIMER SELECTION")
    print("="*80)

    # Helper function for specificity calculation
    def calc_specificity(fseq, rseq, seqs):
        """Specificity = % of non-target sequences NOT matched by primer combo."""
        if not seqs: return None
        ff = re.compile(p2re(fseq), re.I); fr = re.compile(p2re(rc(fseq)), re.I)
        rf = re.compile(p2re(rseq), re.I); rr = re.compile(p2re(rc(rseq)), re.I)
        matches = sum(1 for s in seqs if (ff.search(s) or fr.search(s)) and (rf.search(s) or rr.search(s)))
        return round(100 * (1 - matches / len(seqs)), 1)

    # Full-length sequencing
    print("\n--- Full-Length Sequencing (500-1200 bp amplicons) ---")
    seq_result = select_sequencing(primer_pool, vn_seqs)
    if seq_result:
        a, b, overlap, span, s1, s2, limit = seq_result
        relaxed = f" (overlap relaxed to {limit} bp)" if limit > 200 else ""
        print(f"Span: {a['f']['pos_start']}-{b['r']['pos_end']} ({span} bp), overlap: {overlap} bp{relaxed}")
        print(f"  Combo 1: {a['f']['primer']} + {a['r']['primer']} | {a['amp']} bp | True PCR: {s1}%")
        print(f"  Combo 2: {b['f']['primer']} + {b['r']['primer']} | {b['amp']} bp | True PCR: {s2}%")
        if spec_seqs:
            spec1 = calc_specificity(a['f']['primer'], a['r']['primer'], spec_seqs)
            spec2 = calc_specificity(b['f']['primer'], b['r']['primer'], spec_seqs)
            print(f"  Combo 1 Specificity: {spec1}% | Combo 2 Specificity: {spec2}%")
    else:
        print("No valid overlapping combo found.")

    # qPCR
    print("\n--- qPCR (80-350 bp amplicons) ---")
    qpcr_result = select_qpcr(primer_pool, vn_seqs)
    if qpcr_result:
        c, sens = qpcr_result
        print(f"Best pair: {c['f']['primer']} + {c['r']['primer']}")
        print(f"  Amplicon: {c['amp']} bp | True PCR: {sens}%")
        if spec_seqs:
            spec = calc_specificity(c['f']['primer'], c['r']['primer'], spec_seqs)
            print(f"  Specificity: {spec}%")
    else:
        print("No valid qPCR pair found.")

    # ── Write output ──────────────────────────────────────────────────────────
    lines = [f"# {GENE_NAME} Primer Design & Selection Report\n"]
    lines.append(f"**Design panel:** {total} sequences | **Alignment:** {aln_len} bp | **Min Tm:** {MIN_TM}C\n")
    
    lines.append("## Primer Pool\n")
    header = "| Location | Primer | Tm | GC% | VN% |"
    sep = "|---|---|---|---|---|"
    if cross_seqs:
        header += " " + " | ".join(f"{k}%" for k in cross_seqs.keys()) + " |"
        sep += "|".join(["---"]*len(cross_seqs)) + "|"
    lines.append(header)
    lines.append(sep)
    for p in sorted(primer_pool, key=lambda x: x['pos_start']):
        row = f"| {p['pos_start']}-{p['pos_end']} | `{p['primer']}` | {p['tm']} | {p['gc']*100:.0f} | {p['vn_pct']:.1f} |"
        if cross_seqs:
            row += " " + " | ".join(f"{p['cross_pct'][k]:.1f}" for k in cross_seqs.keys()) + " |"
        lines.append(row)

    lines.append("\n## Selected: Full-Length Sequencing\n")
    if seq_result:
        a, b, overlap, span, s1, s2, limit = seq_result
        lines.append(f"Span: {a['f']['pos_start']}-{b['r']['pos_end']} ({span} bp) | Overlap: {overlap} bp\n")
        lines.append("| Role | Primer | Location | Tm | RC (for PCR) |")
        lines.append("|---|---|---|---|---|")
        lines.append(f"| C1 Fwd | `{a['f']['primer']}` | {a['f']['pos_start']}-{a['f']['pos_end']} | {a['f']['tm']} | (same) |")
        lines.append(f"| C1 Rev | `{a['r']['primer']}` | {a['r']['pos_start']}-{a['r']['pos_end']} | {a['r']['tm']} | `{rc(a['r']['primer'])}` |")
        lines.append(f"| C2 Fwd | `{b['f']['primer']}` | {b['f']['pos_start']}-{b['f']['pos_end']} | {b['f']['tm']} | (same) |")
        lines.append(f"| C2 Rev | `{b['r']['primer']}` | {b['r']['pos_start']}-{b['r']['pos_end']} | {b['r']['tm']} | `{rc(b['r']['primer'])}` |")
        lines.append(f"\n| Amplicon | Size | True PCR VN% |" + (" Specificity% |" if spec_seqs else ""))
        lines.append("|---|---|---|" + ("---|" if spec_seqs else ""))
        if spec_seqs:
            spec1 = calc_specificity(a['f']['primer'], a['r']['primer'], spec_seqs)
            spec2 = calc_specificity(b['f']['primer'], b['r']['primer'], spec_seqs)
            lines.append(f"| C1 | {a['amp']} bp | {s1}% | {spec1}% |")
            lines.append(f"| C2 | {b['amp']} bp | {s2}% | {spec2}% |")
        else:
            lines.append(f"| C1 | {a['amp']} bp | {s1}% |")
            lines.append(f"| C2 | {b['amp']} bp | {s2}% |")
    else:
        lines.append("No valid combo found.\n")

    lines.append("\n## Selected: qPCR\n")
    if qpcr_result:
        c, sens = qpcr_result
        lines.append("| Role | Primer | Location | Tm | RC (for PCR) |")
        lines.append("|---|---|---|---|---|")
        lines.append(f"| Fwd | `{c['f']['primer']}` | {c['f']['pos_start']}-{c['f']['pos_end']} | {c['f']['tm']} | (same) |")
        lines.append(f"| Rev | `{c['r']['primer']}` | {c['r']['pos_start']}-{c['r']['pos_end']} | {c['r']['tm']} | `{rc(c['r']['primer'])}` |")
        spec_str = ""
        if spec_seqs:
            spec = calc_specificity(c['f']['primer'], c['r']['primer'], spec_seqs)
            spec_str = f" | Specificity: {spec}%"
        lines.append(f"\nAmplicon: {c['amp']} bp | True PCR VN%: {sens}%{spec_str}")
    else:
        lines.append("No valid qPCR pair found.\n")

    with open(OUT_FILE, 'w') as fh:
        fh.write('\n'.join(lines))
    print(f"\nOutput saved: {OUT_FILE}")
