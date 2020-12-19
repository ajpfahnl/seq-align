"""
The easiest way to use this package is to `from align import *`.
You will get two classes: 
 * `AlignmentDP` which implements the dynamic programming and 
    backtracking pointers alignment algorithm
 * `AlignmentDC` which implements the divide and conquer
    middle node algorithm for alignment
Useful functions for creating a scoring matrix that both alignment
classes accept are also provided and implemented in `helper.py`.

In general the steps for alignment are:
 0) import the functions

        >>> from align import *

 1) create correctly formatted score matrix:
    
        >>> sc = create_score_matrix_from_file("score_matrix.txt")

        OR

        >>> sc = create_score_matrix_simple(match=5, mismatch=-1, indel=-4)

 2) create Alignment object:

        >>> a = AlignmentDP("AGC", "GCT", sc, local=False)

        OR

        >>> a = AlignmentDC("AGC", "GCT", sc, local=False)

 3) align

        >>> a.align()

 4) print the alignments

        >>> a.print_alignments()

The package also comes with a testing module which can be run with

       $ python3 -m align.test
"""
from .alignmentDP import AlignmentDP
from .alignmentDC import AlignmentDC
from .helper import *