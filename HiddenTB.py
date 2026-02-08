

!pip install biopython

from Bio import SeqIO, pairwise2
from Bio.Seq import Seq
import matplotlib.pyplot as plt

# ---------------------------------------------------------
# Load genome
# ---------------------------------------------------------

def load_genome(fasta_path):
    with open(fasta_path, "r") as f:
        for record in SeqIO.parse(f, "fasta"):
            return str(record.seq).upper()


# ---------------------------------------------------------
# Extract gene using reference coordinates (H37Rv coordinates)
# ---------------------------------------------------------

def extract_gene(genome, start, end):
    return genome[start:end]


# ---------------------------------------------------------
# Align sequences and detect SNPs
# ---------------------------------------------------------

def find_mutations_aligned(ref_seq, query_seq, start_pos):
    alignment = pairwise2.align.globalms(
        ref_seq, query_seq,
        2, -1, -5, -0.5,
        one_alignment_only=True
    )[0]

    aligned_ref, aligned_query = alignment.seqA, alignment.seqB

    mutations = []
    ref_index = start_pos

    for r, q in zip(aligned_ref, aligned_query):
        if r != "-":
            if q != "-" and r != q:
                mutations.append((ref_index, r, q))
            ref_index += 1

    return mutations, len(ref_seq)


# ---------------------------------------------------------
# Paths
# ---------------------------------------------------------

DS_path = r"D:\TB_PROJECT\drug_sensitive sequence\drug_sensitive sequence.fasta"   # H37Rv
DR_path = r"D:\TB_PROJECT\drug_resistant sequence\drug_resistant sequence.fasta"   # SCAID 320.0

DS_genome = load_genome(DS_path)
DR_genome = load_genome(DR_path)

print("Drug-sensitive genome length:", len(DS_genome))
print("Drug-resistant genome length:", len(DR_genome))


# ---------------------------------------------------------
# Resistance genes (H37Rv coordinates)
# ---------------------------------------------------------
resistance_genes = {
    "rpoB":  (759807, 763325),
    "katG":  (2153888, 2156111),
    "inhA":  (2161154, 2162501),
    "gyrA":  (2714915, 2716432),
    "embB":  (4246515, 4249806)
}


# ---------------------------------------------------------
# Analyze genes
# ---------------------------------------------------------
print("\n--- Mutation Analysis (Aligned) ---\n")

mutation_percentages = {}

for gene, (start, end) in resistance_genes.items():
    ref_seq = extract_gene(DS_genome, start, end)
    query_seq = extract_gene(DR_genome, start, end)

    mutations, gene_length = find_mutations_aligned(ref_seq, query_seq, start)

    percent = (len(mutations) / gene_length) * 100
    mutation_percentages[gene] = percent

    print(f"Gene: {gene}")
    print(f"Gene length: {gene_length} bp")
    print(f"Total mutations detected: {len(mutations)}")
    print(f"Mutation percentage: {percent:.2f}%")

    if mutations:
        for pos, ref, alt in mutations[:10]:
            print(f"  {pos}: {ref} → {alt}")
        if len(mutations) > 10:
            print("  ...")
    else:
        print("  No mutations detected")

    print("-" * 40)


# ---------------------------------------------------------
# Visualization: Per-gene mutation density bar chart
# --------------------------------------------------------

genes = list(mutation_percentages.keys())
values = list(mutation_percentages.values())

plt.figure(figsize=(10, 6))
plt.bar(genes, values, color="steelblue", edgecolor="black")

plt.title("Mutation Density per Resistance Gene", fontsize=10)
plt.xlabel("Gene", fontsize=7)
plt.ylabel("Mutation Percentage (%)", fontsize=10)

plt.grid(axis="y", linestyle="--", alpha=0.5)
plt.tight_layout()
plt.show()
