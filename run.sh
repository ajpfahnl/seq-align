#!/bin/bash

echo "This program was tested with Python 3.8.6."

# check for required python modules
(python3 -m memory_profiler --version && python3 -c 'import numpy') >/dev/null 2>&1 || 
{ echo '
WARNING: Either numpy or memory-profiler are not installed.
         Please run setup.sh to install necessary packages.
'; exit 1; }

echo
echo 'For documentation on the align package, please read through `__init__.py` and the docstrings in the alignment classes.'
echo
echo ">> A Simple Example: "

align_script="
from align import *

print('Dynamic programming: global')
sc = create_score_matrix_from_file('score_matrix.txt')
a = AlignmentDP('AGC', 'GCT', sc, local=False)
a.align()
a.print_alignments()

print('Divide and conquer: local')
sc = create_score_matrix_simple(match=1, mismatch=-1, indel=-1)
b = AlignmentDC('AGC', 'GCT', sc, local=True)
b.align()
b.print_alignments()
"

echo "--- Python Script ---"
echo "$align_script"
echo "------ Output -------"
echo

echo "$align_script" | python3

echo
read -p "Press ENTER to proceed with testing..."

# run test module
python3 -m align.test
