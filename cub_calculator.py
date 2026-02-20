import os
import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from collections import defaultdict
import numpy as np
import concurrent.futures
import warnings
warnings.filterwarnings("ignore")

# === PARAMETERS ===
input_dir = "ffn_files"  # Folder containing multifasta gene files
output_dir = "genomes_csv"
os.makedirs(output_dir, exist_ok=True)

NUM_WORKERS = 30  # Number of threads for parallel processing

# Standard Genetic Code Dictionary
codon_table = {
    'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L',
    'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L',
    'ATT':'I', 'ATC':'I', 'ATA':'I', 'ATG':'M',
    'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V',
    'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S',
    'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P',
    'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T',
    'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
    'TAT':'Y', 'TAC':'Y', 'TAA':'*', 'TAG':'*',
    'CAT':'H', 'CAC':'H', 'CAA':'Q', 'CAG':'Q',
    'AAT':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K',
    'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E',
    'TGT':'C', 'TGC':'C', 'TGA':'*', 'TGG':'W',
    'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R',
    'AGT':'S', 'AGC':'S', 'AGA':'R', 'AGG':'R',
    'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G',
}

def calculate_rscu(codon_counts):
    """Calculates Relative Synonymous Codon Usage (RSCU)."""
    rscu = {}
    aa_to_codons = defaultdict(list)
    for codon, aa in codon_table.items():
        aa_to_codons[aa].append(codon)
    for aa, codons in aa_to_codons.items():
        total = sum([codon_counts.get(c, 0) for c in codons])
        for codon in codons:
            observed = codon_counts.get(codon, 0)
            expected = total / len(codons) if len(codons) > 0 else 0
            rscu[codon] = observed / expected if expected > 0 else 0
    return rscu

def calculate_enc(seq):
    """Calculates Effective Number of Codons (ENC)."""
    codon_counts = defaultdict(int)
    if len(seq) % 3 != 0:
        return np.nan
    for i in range(0, len(seq), 3):
        codon = seq[i:i+3]
        codon_counts[codon] += 1
    N = sum(codon_counts.values())
    if N == 0:
        return np.nan
    freqs = np.array(list(codon_counts.values())) / N
    sum_sq_freqs = sum(freqs ** 2)
    return 1 / sum_sq_freqs if sum_sq_freqs > 0 else np.nan

def calculate_codon_usage_profile(seqs):
    """Generates a codon usage profile from a set of reference sequences."""
    codon_freqs = defaultdict(int)
    total_codons = 0
    for seq in seqs:
        seq = str(seq).upper()
        for i in range(0, len(seq) - 2, 3):
            codon = seq[i:i+3]
            if codon in codon_table:
                codon_freqs[codon] += 1
                total_codons += 1
    profile = {}
    for codon in codon_table:
        profile[codon] = codon_freqs[codon] / total_codons if total_codons > 0 else 0
    return profile

def calculate_cai(seq, ref_profile):
    """Calculates Codon Adaptation Index (CAI) relative to a reference profile."""
    seq = seq.upper()
    weights = []
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]
        if codon in codon_table and ref_profile.get(codon, 0) > 0:
            aa = codon_table[codon]
            synonymous_codons = [c for c, a in codon_table.items() if a == aa]
            max_freq = max([ref_profile.get(c, 0) for c in synonymous_codons])
            if max_freq > 0:
                w = ref_profile.get(codon, 0) / max_freq
                weights.append(w)
    if len(weights) == 0:
        return np.nan
    return np.exp(np.mean(np.log(weights)))

def calculate_genome_rscu_profile(records):
    """Calculates the global RSCU profile for the entire genome."""
    codon_counts_total = defaultdict(int)
    total_codons = 0
    for r in records:
        seq = str(r.seq).upper()
        for i in range(0, len(seq)-2, 3):
            codon = seq[i:i+3]
            if codon in codon_table:
                codon_counts_total[codon] += 1
                total_codons += 1
    genome_profile = {codon: codon_counts_total[codon] / total_codons if total_codons > 0 else 0 for codon in codon_table}
    return genome_profile

def rscu_distance(rscu_gene, genome_rscu_profile):
    """Calculates Mean Absolute Difference between gene RSCU and genome profile."""
    diffs = []
    for codon in codon_table:
        g_val = genome_rscu_profile.get(codon, 0)
        gene_val = rscu_gene.get(codon, 0)
        diffs.append(abs(gene_val - g_val))
    return np.mean(diffs)

def process_genome(filename):
    """Processes a single genome file and extracts CUB metrics for all genes."""
    genome_id = filename.replace(".ffn", "")
    filepath = os.path.join(input_dir, filename)
    records = list(SeqIO.parse(filepath, "fasta"))

    # Select highly expressed reference genes (Ribosomal proteins, etc.)
    ref_seqs = [r.seq for r in records if any(k in r.description.lower() for k in 
                ["ribosomal", "50s", "30s", "ssra", "rimm", "rimp", "rna small subunit"])]
    
    ref_profile = calculate_codon_usage_profile(ref_seqs)
    genome_rscu_profile = calculate_genome_rscu_profile(records)

    gene_data = []
    gc_contents = []
    rscu_dists = []

    for record in records:
        desc = record.description.lower()

        # Handle non-coding RNA or other exclusions
        if any(tag in desc for tag in ["trna", "rrna", "non-coding", "ncrna", "misc_rna"]):
            gc_val = gc_fraction(record.seq) * 100
            gene_data.append({
                "genome": genome_id,
                "gene_id": record.id,
                "description": record.description,
                "ENC": np.nan,
                "CAI": np.nan,
                "GC_content": gc_val,
                "GC_zscore_self": np.nan,
                "Codon_Avg_Dist_self": np.nan,
                "Quality_note": "Excluded: non-coding RNA"
            })
            continue

        seq = str(record.seq).upper()
        if len(seq) < 30:
            continue

        quality_note = "OK"
        if len(seq) % 3 != 0:
            quality_note = "Sequence not valid (not multiple of 3)"
        elif len(seq) // 3 < 100:
            quality_note = "Short gene (<100 codons)"

        if len(seq) % 3 == 0:
            codon_counts = defaultdict(int)
            for i in range(0, len(seq), 3):
                codon = seq[i:i+3]
                codon_counts[codon] += 1
            rscu = calculate_rscu(codon_counts)
            enc = calculate_enc(seq)
            cai = calculate_cai(seq, ref_profile)
        else:
            rscu = {codon: np.nan for codon in codon_table}
            enc = np.nan
            cai = np.nan

        gc_val = gc_fraction(seq) * 100
        dist_val = rscu_distance(rscu, genome_rscu_profile)

        gc_contents.append(gc_val)
        rscu_dists.append(dist_val)

        gene_data.append({
            "genome": genome_id,
            "gene_id": record.id,
            "description": record.description,
            "ENC": enc,
            "CAI": cai,
            "GC_content": gc_val,
            "GC_zscore_self": np.nan,
            "Codon_Avg_Dist_self": dist_val,
            "Quality_note": quality_note,
            **rscu
        })

    # Calculate GC Z-scores based on the genome's own distribution
    if gc_contents:
        gc_mean = np.mean(gc_contents)
        gc_std = np.std(gc_contents)
        for gene in gene_data:
            if pd.notna(gene["GC_content"]):
                gene["GC_zscore_self"] = (gene["GC_content"] - gc_mean) / gc_std if gc_std > 0 else 0

    df = pd.DataFrame(gene_data)
    df.to_csv(os.path.join(output_dir, f"{genome_id}_CUB.csv"), index=False)
    return gene_data

# === MAIN EXECUTION ===
if __name__ == "__main__":
    all_genes_data = []

    ffn_files = [f for f in os.listdir(input_dir) if f.endswith(".ffn")]
    total_files = len(ffn_files)
    completed = 0

    print(f"Starting analysis of {total_files} genomes using {NUM_WORKERS} workers...")

    with concurrent.futures.ProcessPoolExecutor(max_workers=NUM_WORKERS) as executor:
        futures = {executor.submit(process_genome, f): f for f in ffn_files}
        for future in concurrent.futures.as_completed(futures):
            result = future.result()
            all_genes_data.extend(result)
            completed += 1
            print(f"Progress: {completed}/{total_files} files completed")

    print("\n Processing finished. All CSV files are located in:", output_dir)
