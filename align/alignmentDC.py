import numpy as np
import math
from .alignmentDP import AlignmentDP

class AlignmentDC():
    """
    Divide and conquer alignment class.
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
        created by align()
    local: bool
        whether to locally align or not

    Methods
    -------
    align():
        calls on align_helper() to recursively align locally or globally based on the value of local
        returns alignments found
    print_alignments():
        print the alignments in a pretty way
    """
    def __init__(self, s1, s2, score_matrix, local=False):
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
        """
        self.s1 = s1
        self.s2 = s2
        self.score_matrix = score_matrix
        self.local = local

        # nucleotide_index
        self.ni = {"A": 0, "C": 1, "G": 2, "T": 3, "-": 4}

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

    def score(self, X,Y, max_bool=False, max_val=0):
        """
        Returns the last column in the alignment scoring matrix between X and Y
        for global alignment. Also can return where a certain value occurred.
        """
        if max_bool:
            max_locs = []
        s0 = np.zeros(len(Y)+1)
        s1 = np.zeros(len(Y)+1)
        for j in range(1, len(Y)+1):
            s0[j] = s0[j-1] + self.indel(Y[j-1])
        for i in range(1, len(X)+1):
            s1[0] = s0[0] + self.indel(X[i-1])
            for j in range(1, len(Y)+1):
                score_match = s0[j-1] + self.match(X[i-1], Y[j-1])
                score_del = s0[j] + self.indel(X[i-1])
                score_ins = s1[j-1] + self.indel(Y[j-1])
                s1[j] = max(score_match, score_del, score_ins)
                if max_bool and s1[j] == max_val:
                        max_locs.append((i,j))
            s0 = np.copy(s1)
        if max_bool:
            return s1, max_locs, max_val
        return s1
    
    def score_local(self, X,Y, max_bool=False, max_val=0, soft_max=True):
        """
        Returns the last column in the alignment scoring matrix between X and Y
        for local alignment. Also can return the maximum value and where it occurred.
        """
        if max_bool:
            max_locs = []
        s0 = np.zeros(len(Y)+1)
        s1 = np.zeros(len(Y)+1)
        for i in range(1, len(X)+1):
            s1[0] = 0
            for j in range(1, len(Y)+1):
                score_match = s0[j-1] + self.match(X[i-1], Y[j-1])
                score_del = s0[j] + self.indel(X[i-1])
                score_ins = s1[j-1] + self.indel(Y[j-1])
                s1[j] = max(0, score_match, score_del, score_ins)
                if max_bool:
                    if soft_max and s1[j] > max_val:
                        max_val = s1[j]
                        max_locs = [(i, j)]
                    elif s1[j] == max_val:
                        max_locs.append((i,j))
            s0 = np.copy(s1)
        if max_bool:
            return s1, max_locs, max_val
        return s1
    
    def align_helper(self, X, Y, local=False):
        """
        Divide and conquer alignment algorithm

        Parameters
        ----------
        X : string
            first string to align
        Y : string
            second string to align
        local : bool, optional
            set to True for local alignment
            defaults to global alignment

        Returns
        -------
        list
            a list of tuples with the string representations of the alignments between X and Y
        """
        Z = ""
        W = ""
        xlen = len(X)
        ylen = len(Y)

        # Base cases:
        ## case 1: X is ""
        if xlen == 0:
            for i in range(ylen):
                Z += '-'
                W += Y[i]
            ZWs = [(Z, W)]
        ## case 2: Y is ""
        elif ylen == 0:
            for i in range(xlen):
                Z = Z + X[i]
                W = W + '-'
            ZWs = [(Z, W)]
        ## case 3: one of X or Y is a single character
        elif xlen == 1 or ylen == 1:
            ZWs = AlignmentDP(X, Y, self.score_matrix, local=local).align()
        
        # Recursive alignment:
        else:
            xmid = len(X)//2

            # find middle node(s)
            if local:
                scoreL = self.score_local(X[:xmid], Y)
                scoreR = self.score_local(X[xmid:][::-1], Y[::-1])
            else:
                scoreL = self.score(X[:xmid], Y)
                scoreR = self.score(X[xmid:][::-1], Y[::-1])
            score = scoreL + np.flip(scoreR)
            ymids = np.where(score == np.max(score))[0]

            del score
            del scoreR
            del scoreL

            # find alignment(s)
            ZWs = set()
            for ymid in ymids:
                ZWLs = self.align_helper(X[:xmid], Y[:ymid])
                ZWRs = self.align_helper(X[xmid:], Y[ymid:])
                for ZWL in ZWLs:
                    # ZL = ZWL[0]
                    # WL = ZWL[1]
                    for ZWR in ZWRs:
                        # ZR = ZWR[0]
                        # WR = ZWR[1]
                        ZWs.add((ZWL[0]+ZWR[0], ZWL[1]+ZWR[1]))
        return list(ZWs)

    def align(self):
        """
        Run the recursive alignment algorithm and return distinct alignments

        Returns
        -------
        list
            a list of tuples with the string representations of the alignments between s1 and s2
        """
        #
        # global
        #
        if not self.local:
            alignments = self.align_helper(self.s1, self.s2)
            alignments = list(set(alignments))
            alignments.sort()
            self.alignments = alignments
            return alignments
        
        #
        # local
        #
        alignments = []
        # find the alignment end nodes
        _, max_loc_ends, maxval = self.score_local(self.s1, self.s2, max_bool=True)
        # find the alignment start nodes associated with each end node
        max_loc_starts = []
        for mle in max_loc_ends:
            _, mlss, _ = self.score(self.s1[:mle[0]][::-1], self.s2[:mle[1]][::-1], max_bool=True, max_val=maxval)
            max_loc_starts.append([(mle[0]-mls[0], mle[1]-mls[1]) for mls in mlss])
        # divide and conquer on each start node and end node pair and add the alignments found
        for max_loc_end, mlss in zip(max_loc_ends, max_loc_starts):
            for max_loc_start in mlss:
                s1_trim = self.s1[max_loc_start[0]:max_loc_end[0]]
                s2_trim = self.s2[max_loc_start[1]:max_loc_end[1]]
                alignments.extend(self.align_helper(s1_trim, s2_trim))
        alignments.sort()
        self.alignments = alignments
        return alignments

    def print_alignments(self):
        """
        Prints the optimal alignments found by align
        """
        for i, alignment in enumerate(self.alignments):
            print(f"   alignment {i+1}:")
            print(f"   {alignment[0]}")
            print(f"   {alignment[1]}")
        print(f"   {len(self.alignments)} alignment(s) total.")