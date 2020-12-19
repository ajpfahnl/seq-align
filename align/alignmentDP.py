import numpy as np

class AlignmentDP():
    """
    Dynamic programming alignment class.
    Attributes
    ----------
    s1 : str
        first string to align
    s2 : str
        second string to align
    score_matrix : numpy array
        indel, mismatch/match scores
    alignments : list of tuples
        the elments of the tuples are string representations of the alignments of s1 and s2, respectively
        created by align(), specifically find_alignments()
    local: bool
        whether to locally align or not

    Methods
    -------
    align():
        runs create_matrix_and_backpointers() and find_alignments()
        returns alignments found
    create_matrix_and_backpointers():
        creates alignment matrix and backpointers
    find_alignments():
        finds alignments globally or locally depending on the value of local
        requires that create_matrix_and_backpointers() is run first
    print_alignment_matrix():
        print the alignment matrix
    print_backtrack_matrix();
        print the backtracking matrix in a pretty way
    print_alignments():
        print the alignments in a pretty way
    """
    def __init__(self, s1, s2, score_matrix, local=False, verbose=False):
        """
        Parameters
        ----------
        s1 : string
            first string to align
        s2 : string
            second string to align
        score_matrix :
            numpy array with shape (5, 5)
            indices 0 to 4 for rows and columns follow A, C, G, T, -
            should be a symmetric matrix
        local : bool, optional
            if True, changes to local alignment
            default is global alignment
        verbose : bool, optional
            step through the alignment and backtracking pointer creation process
        """
        self.s1 = "-" + s1
        self.s2 = "-" + s2
        self.score_matrix = score_matrix
        self.local = local
        self.verbose = verbose

        # nucleotide_index
        self.ni = {"A": 0, "C": 1, "G": 2, "T": 3, "-": 4}

        # backtracking pointer index
        self.bi = {"→": 0, "↓": 1, "↘": 2}

    def indel(self, C):
        """
        Indel score (usually a negative penalty)
        """
        return self.score_matrix[self.ni[C], self.ni['-']]
    def match(self, C1, C2):
        """
        Match/mismatch score (usually a negative penalty if mismatch)
        """
        return self.score_matrix[self.ni[C1], self.ni[C2]]

    def create_matrix_and_backpointers(self):
        """
        Create alignment matrix and backpointers

        Returns
        -------
        (numpy matrix, numpy matrix):
            Returns the alignment matrix and the backpointers matrix
        """
        verbose = self.verbose
        local = self.local
        s1len = len(self.s1)
        s2len = len(self.s2)
        
        # allocate alignment matrix and backpointer matrix
        m = np.array([[0.0 for i in range(s1len)] for j in range(s2len)])
        b = np.array([[(False, False, False) for i in range(s1len)] for j in range(s2len)])

        # initialization for global
        if not local:
            # initialize first row
            for j in range(1,s1len):
                m[0, j] = m[0, j-1] + self.indel(self.s1[j])
                b[0, j, self.bi["→"]] = True

            # initialize first column
            for i in range(1,s2len):
                m[i, 0] = m[i-1, 0] + self.indel(self.s2[i])
                b[i, 0, self.bi["↓"]] = True
            
        # fill the rest of the alignment scores
        for i in range(1, s2len):
            for j in range(1, s1len):
                
                match_mismatch = m[i-1, j-1] + self.match(self.s2[i], self.s1[j])
                indel_s1_gap = m[i-1, j] + self.indel(self.s2[i])
                indel_s2_gap = m[i, j-1] + self.indel(self.s1[j])

                if local:
                    m[i, j] = max(0, match_mismatch, indel_s1_gap, indel_s2_gap)
                else:
                    m[i, j] = max(match_mismatch, indel_s1_gap, indel_s2_gap)

                if m[i, j] == match_mismatch: 
                    b[i, j, self.bi["↘"]] = True
                if m[i, j] == indel_s1_gap:
                    b[i, j, self.bi["↓"]] = True
                if m[i, j] == indel_s2_gap:
                    b[i, j, self.bi["→"]] = True

                if verbose:
                    print(f"{self.s1[j]} {self.s2[i]}")
                    print(m)
                    print("Press ENTER to continue...")
                    input()
        self.m = m
        self.b = b
        return m, b

    def print_alignment_matrix(self):
        """
        Print the alignment matrix. 
        
        Run create_matrix_and_backpointers first.
        """
        print(self.m)

    def print_backtrack_matrix(self, length=3):
        """
        Print the backtracking pointers matrix. 
        
        Run create_matrix_and_backpointers first.
        """
        for row in self.b:
            for pointers in row:
                to_print = ''
                if pointers[0] == True:
                    to_print += "→"
                if pointers[1] == True:
                    to_print += "↓"
                if pointers[2] == True:
                    to_print += "↘"
                if not np.any(pointers):
                    to_print += "s"
                print('{:{length}}'.format(to_print, length=length), end='')
            print()

    def find_alignments(self):
        """
        Find all optimal alignments.

        Run create_matrix_and_backpointers first.

        Returns
        -------
        list
            a list of tuples with the string representations of the alignments between s1 and s2
        """
        local = self.local
        # initalize alignment stack: this will hold intermediate alignments until
        # a full alignment is found and placed in final_alignment
        i, j = self.m.shape
        if local:
            # find largest value (can have multiple) and backtrack from there
            indices = np.where(self.m==np.max(self.m))
            alignment_stack = []
            for i in range(len(indices[0])):
                alignment_stack.append(((indices[0][i],indices[1][i]),'', ''))

        else:
            alignment_stack = [((i-1,j-1),'', '')]
        alignments = []

        def process_new_alignment(new_location, new_s2_alignment, new_s1_alignment):
            if not np.any(self.b[new_location[0], new_location[1]]):
                alignments.append((new_s1_alignment[::-1], new_s2_alignment[::-1]))
            else:
                alignment_stack.append((new_location, new_s2_alignment, new_s1_alignment))
        
        # find alignments
        while len(alignment_stack) != 0:
            location, s2_alignment, s1_alignment = alignment_stack.pop()

            # match/mismatch
            if self.b[location[0], location[1], 2]:
                new_s2_alignment = s2_alignment + self.s2[location[0]]
                new_s1_alignment = s1_alignment + self.s1[location[1]]
                new_location = (location[0]-1, location[1]-1)
                process_new_alignment(new_location, new_s2_alignment, new_s1_alignment)
            # gap in s1
            if self.b[location[0], location[1], 1]:
                new_s2_alignment = s2_alignment + self.s2[location[0]]
                new_s1_alignment = s1_alignment + '-'
                new_location = (location[0]-1, location[1])
                process_new_alignment(new_location, new_s2_alignment, new_s1_alignment)
            # gap in s2
            if self.b[location[0], location[1], 0]:
                new_s2_alignment = s2_alignment + '-'
                new_s1_alignment = s1_alignment + self.s1[location[1]]
                new_location = (location[0], location[1]-1)
                process_new_alignment(new_location, new_s2_alignment, new_s1_alignment)
        alignments.sort()
        self.alignments = alignments
        return alignments

    def print_alignments(self):
        """
        Prints the optimal alignments found by find_alignments.
        """
        for i, alignment in enumerate(self.alignments):
            print(f"   alignment {i+1}:")
            print(f"   {alignment[0]}")
            print(f"   {alignment[1]}")
        print(f"   {len(self.alignments)} alignment(s) total.")

    def align(self):
        """
        Align s1 and s2

        Returns
        -------
        list
            a list of tuples with the string representations of the alignments between s1 and s2
        """
        self.create_matrix_and_backpointers()
        return self.find_alignments()
