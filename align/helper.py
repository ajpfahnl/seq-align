import numpy as np
def create_score_matrix_simple(match=1, mismatch=-1, indel=-1):
    """
    Creates a score matrix from three scores for matches, mismatches, and indels.

    Parameters
    ----------
    match : int, default 1
        score for matches (usually positive)
    mismatch : int, default -1
        score for mismatches (usually negative)
    indel : int, default -1
        score for indels (usually negative)

    Returns
    -------
    numpy array with shape (5, 5)
        returns a symmetric matrix
        indices 0 to 4 for rows and columns follow A, C, G, T, -
    """
    score_matrix = [[match, mismatch, mismatch, mismatch, indel],
            [mismatch, match, mismatch, mismatch, indel],
            [mismatch, mismatch, match, mismatch, indel],
            [mismatch, mismatch, mismatch, match, indel],
            [indel, indel, indel, indel, 0.0],
            ]
    return np.array(score_matrix)

def create_score_matrix_from_file(file_path):
    """
    Creates a score matrix from a file.

    File contents must have the same format as the example below:

            A       C       G       T       -
    A       1       -1      -2      -3      -4
    C       -1      6       -3      -4      -5
    G       -2      -3      5       -5      -6
    T       -3      -4      -5      4       -1
    -       -4      -5      -6      -1

    Notes: 
     * (-, -) does not have a score
     * all values are separated by some form of whitespace
     * no empty rows before the start of the matrix

    Parameters
    ----------
    file_path : string
        path to the score matrix file
        can be a relative path from where you are executing this function OR
        an absolute path

    Returns
    -------
    numpy array with shape (5, 5)
        returns a symmetric matrix
        indices 0 to 4 for rows and columns follow A, C, G, T, -
    """
    score_matrix =[]
    with open(file_path, 'r') as f:
        f.readline()
        for i, line in enumerate(f):
            if i == 5:
                break
            line = line.strip().split()
            score_matrix.append([int(s) for s in line[1:]])
    score_matrix[-1].append(0.0)
    return np.array(score_matrix)
