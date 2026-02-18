# Importing Dependencies
import pytest
from align import NeedlemanWunsch, read_fasta
import numpy as np

def test_nw_alignment():
    """
    TODO: Write your unit test for NW alignment
    using test_seq1.fa and test_seq2.fa by
    asserting that you have correctly filled out
    the your 3 alignment matrices.
    Use the BLOSUM62 matrix and a gap open penalty
    of -10 and a gap extension penalty of -1.
    """
    seq1, _ = read_fasta("./data/test_seq1.fa")  # MYQR
    seq2, _ = read_fasta("./data/test_seq2.fa")  # MQR

    nw = NeedlemanWunsch(
        sub_matrix_file="./substitution_matrices/BLOSUM62.mat",
        gap_open=-10,
        gap_extend=-1
    )
    nw.align(seq1, seq2)

    # Check some key values in the align matrix (M)
    # M[0,0] should be 0 (base case)
    assert nw._align_matrix[0, 0] == 0

    # M[1,1] = M-M match = 0 + 5 = 5
    assert nw._align_matrix[1, 1] == 5

    # M[4,3] is the final match score (R-R after optimal path)
    assert nw._align_matrix[4, 3] == 4

    # Check gap matrices initialization
    # First row of gapA should have gap penalties
    assert nw._gapA_matrix[0, 1] == -11  # gap_open + 1*gap_extend
    assert nw._gapA_matrix[0, 2] == -12
    assert nw._gapA_matrix[0, 3] == -13

    # First column of gapB should have gap penalties
    assert nw._gapB_matrix[1, 0] == -11
    assert nw._gapB_matrix[2, 0] == -12
    

def test_nw_backtrace():
    """
    TODO: Write your unit test for NW backtracing
    using test_seq3.fa and test_seq4.fa by
    asserting that the backtrace is correct.
    Use the BLOSUM62 matrix. Use a gap open
    penalty of -10 and a gap extension penalty of -1.
    """
    seq3, _ = read_fasta("./data/test_seq3.fa")  # MAVHQLIRRP
    seq4, _ = read_fasta("./data/test_seq4.fa")  # MQLIRHP

    nw = NeedlemanWunsch(
        sub_matrix_file="./substitution_matrices/BLOSUM62.mat",
        gap_open=-10,
        gap_extend=-1
    )
    score, align1, align2 = nw.align(seq3, seq4)

    # Check alignment score matches expected value from README
    assert score == 17

    # Check the aligned sequences match expected output
    assert align1 == "MAVHQLIRRP"
    assert align2 == "M---QLIRHP"




