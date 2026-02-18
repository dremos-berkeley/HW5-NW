# Import NeedlemanWunsch class and read_fasta function
from align import read_fasta, NeedlemanWunsch

def main():
    """
    This function should
    (1) Align all species to humans and print species in order of most similar to human BRD
    (2) Print all alignment scores between each species BRD2 and human BRD2
    """
    hs_seq, hs_header = read_fasta("./data/Homo_sapiens_BRD2.fa")
    gg_seq, gg_header = read_fasta("./data/Gallus_gallus_BRD2.fa")
    mm_seq, mm_header = read_fasta("./data/Mus_musculus_BRD2.fa")
    br_seq, br_header = read_fasta("./data/Balaeniceps_rex_BRD2.fa")
    tt_seq, tt_header = read_fasta("./data/tursiops_truncatus_BRD2.fa")

    # Create aligner with BLOSUM62 and specified gap penalties
    nw = NeedlemanWunsch(
        sub_matrix_file="./substitution_matrices/BLOSUM62.mat",
        gap_open=-10,
        gap_extend=-1
    )

    # Align each species to human and store results
    species_scores = []

    score, _, _ = nw.align(hs_seq, gg_seq)
    species_scores.append(("Gallus gallus", score))

    score, _, _ = nw.align(hs_seq, mm_seq)
    species_scores.append(("Mus musculus", score))

    score, _, _ = nw.align(hs_seq, br_seq)
    species_scores.append(("Balaeniceps rex", score))

    score, _, _ = nw.align(hs_seq, tt_seq)
    species_scores.append(("Tursiops truncatus", score))

    # Sort by score descending (most similar first)
    species_scores.sort(key=lambda x: x[1], reverse=True)

    # Print species in order of similarity
    print("Species in order of most similar to least similar to human BRD2:")
    for species, _ in species_scores:
        print(f"  {species}")

    # Print alignment scores
    print("\nAlignment scores:")
    for species, score in species_scores:
        print(f"  {species}: {score}")
    

if __name__ == "__main__":
    main()
