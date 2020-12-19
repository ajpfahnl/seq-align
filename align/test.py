import numpy as np
from .alignmentDP import AlignmentDP
from .alignmentDC import AlignmentDC
from .helper import create_score_matrix_simple
from unittest import TestCase
import unittest
import math
import random
import time
from statistics import mean

MP_AVAIL = True
try:
    from memory_profiler import memory_usage
except:
    MP_AVAIL = False

class AlignTest(TestCase):
    """
    Testing base class with alignment functions and assertion checks

    Test cases are named with the "test" prefix.
    """
    def __init__(self, *args, **kwargs):
        super(AlignTest, self).__init__(*args, **kwargs)
        self.print_bool = True

    def create_score_matrix(self, match=1, mismatch=-1, indel=-1):
        return(create_score_matrix_simple(match, mismatch, indel))

    def find_alignments_dp(self, s1, s2, score_matrix, local=False, ret=True):
        if self.print_bool:
            print(">> Dynamic Programming")
        align = AlignmentDP(s1, s2, score_matrix, local=local)
        dp_align = align.align()
        if self.print_bool:
            align.print_alignments()
        if ret == True:
            return dp_align

    def find_alignments_dc(self, s1, s2, score_matrix, local=False, ret=True):
        if self.print_bool:
            print(">> Divide and Conquer")
        align = AlignmentDC(s1, s2, score_matrix, local=local)
        dc_align = align.align()
        if self.print_bool:
            align.print_alignments()
        if ret == True:
            return dc_align

    def compare_dp_dc(self, dp_align, dc_align):
        self.assertEqual(len(dp_align), len(dc_align), msg="Number of alignments found differ. Check output.")
        for alignment_pair in zip(dp_align, dc_align):
            self.assertEqual(alignment_pair[0][0], alignment_pair[1][0], msg="Inconsistent alignment. Check output.")
            self.assertEqual(alignment_pair[0][1], alignment_pair[1][1], msg="Inconsistent alignment. Check output.")
        if self.print_bool:
            print(">> Alignments match")

    def find_alignments(self, s1, s2, score_matrix, local=False):
        if self.print_bool:
            if local:
                print(">> Local alignment")
            else:
                print(">> Global alignment")
        dp_align = self.find_alignments_dp(s1, s2, score_matrix, local)
        dc_align = self.find_alignments_dc(s1, s2, score_matrix, local)
        self.compare_dp_dc(dp_align, dc_align)
    
    def print_originals(self, s1, s2):
        if self.print_bool:
            print(f">> Original strings")
            print(f"   {s1}")
            print(f"   {s2}")

    def print_scores(self, score_matrix):
        if self.print_bool:
            print(f">> score matrix")
            print(score_matrix)



class CorrectnessTest_a_Simple(AlignTest):
    """
    Simple correctness tests
    """
    def check_correct(self, s1, s2, score_matrix, local=False, msg=None):
        if self.print_bool:
            print()
            if msg:
                print(f"{msg}")
        self.print_originals(s1, s2)
        self.print_scores(score_matrix)
        self.find_alignments(s1, s2, score_matrix, local=local)

    def random_score_matrix(self, seed=1, low=1, high=21):
        np.random.seed(seed)
        score_matrix = np.random.randint(-high, -low, (5,5))
        for i in range(5):
            score_matrix[i,i] = np.abs(score_matrix[i,i])
        for i in range(5):
            for j in range(5):
                score_matrix[i,j] = score_matrix[j,i]
        score_matrix[4,4] = 0.0
        return score_matrix
    
    def test_1a_global(self):
        s1 = "CTATGCCA"
        s2 = "CCTACA"
        score_matrix = self.create_score_matrix(match=2, mismatch=-1, indel=-1)
        self.check_correct(s1, s2, score_matrix)

    def test_1a_local(self):
        s1="CTATGCCA"
        s2="CCTACA"
        score_matrix = self.create_score_matrix(match=2, mismatch=-1, indel=-1)
        self.check_correct(s1, s2, score_matrix, local=True)

    def test_1b_global(self):
        s1="CTATGCCA"
        s2="CCTACA"
        score_matrix = self.create_score_matrix(match=2, mismatch=-1, indel=-5)
        self.check_correct(s1, s2, score_matrix)

    def test_1b_local(self):
        s1="CTATGCCA"
        s2="CCTACA"
        score_matrix = self.create_score_matrix(match=2, mismatch=-1, indel=-5)
        self.check_correct(s1, s2, score_matrix, local=True)

    def test_2_global(self):
        s1="CTATGCCA"
        s2="CTATGCCA"
        score_matrix = self.create_score_matrix(match=2, mismatch=-5, indel=-5)
        self.check_correct(s1, s2, score_matrix)

    def test_2_local(self):
        s1="CTATGCCA"
        s2="CTATGCCA"
        score_matrix = self.create_score_matrix(match=2, mismatch=-5, indel=-5)
        self.check_correct(s1, s2, score_matrix, local=True)

    def test_3_global(self):
        s1="AAAAAAAA"
        s2="AAAAAAAA"
        score_matrix = self.create_score_matrix(match=2, mismatch=-5, indel=-5)
        self.check_correct(s1, s2, score_matrix)

    def test_3_local(self):
        s1="AAAAAAAA"
        s2="AAAAAAAA"
        score_matrix = self.create_score_matrix(match=2, mismatch=-5, indel=-5)
        self.check_correct(s1, s2, score_matrix, local=True)

    def test_4a_global(self):
        s1="AAAAGTCAAAA"
        s2="GTC"
        score_matrix = self.create_score_matrix(match=2, mismatch=-5, indel=-5)
        self.check_correct(s1, s2, score_matrix)

    def test_4a_local(self):
        s1="AAAAGTCAAAA"
        s2="GTC"
        score_matrix = self.create_score_matrix(match=2, mismatch=-5, indel=-5)
        self.check_correct(s1, s2, score_matrix, local=True)

    def test_4b_global(self):
        s1="GTC"
        s2="AAAAGTCAAAA"
        score_matrix = self.create_score_matrix(match=2, mismatch=-5, indel=-5)
        self.check_correct(s1, s2, score_matrix)

    def test_4b_local(self):
        s1="GTC"
        s2="AAAAGTCAAAA"
        score_matrix = self.create_score_matrix(match=2, mismatch=-5, indel=-5)
        self.check_correct(s1, s2, score_matrix, local=True)

    def test_5_global(self):
        s1="AAAAGTCAAAAATGAAAAA"
        s2="GTCTGA"
        score_matrix = self.create_score_matrix(match=2, mismatch=-5, indel=-5)
        self.check_correct(s1, s2, score_matrix)

    def test_5_local(self):
        s1="AAAAGTCAAAAATGAAAAA"
        s2="GTCTGA"
        score_matrix = self.create_score_matrix(match=2, mismatch=-5, indel=-5)
        self.check_correct(s1, s2, score_matrix, local=True)

    def test_random_score_global(self):
        if self.print_bool:
            print()
        s1 = "CTATGCCA"
        s2 = "CCTACA"
        self.print_originals(s1, s2)
        for seed in [1,5,10]:
            score_matrix = self.random_score_matrix(seed)
            self.print_scores(score_matrix)
            self.find_alignments(s1, s2, score_matrix, local=False)

    def test_random_score_local(self):
        if self.print_bool:
            print()
        s1 = "CTATGCCA"
        s2 = "CCTACA"
        self.print_originals(s1, s2)
        for seed in [1,5,10]:
            score_matrix = self.random_score_matrix(seed)
            self.print_scores(score_matrix)
            self.find_alignments(s1, s2, score_matrix, local=True)


@unittest.skipIf(not MP_AVAIL, "memory-map module not installed")
class MemTimeTest(AlignTest):
    """
    Memory and time check tests with their helper functions.
    """
    def __init__(self, *args, **kwargs):
        super(MemTimeTest, self).__init__(*args, **kwargs)
        self.letters = ('A', 'C', 'G', 'T')
        self.print_bool = False
    
    def random_string(self, length, seed=10):
        random.seed(seed)
        s = ""
        for i in range(length):
            s += random.choice(self.letters)
        return s

    def time_mem_sanity(self, s1, s2, score_matrix, length=None, local=False, print_bool=True):
        length_sum = len(s1) + len(s2)
        if length == None:
            length = length_sum//2
        # time and memory check
        tic_dc = time.perf_counter()
        mem_usage_dc, dc_align = memory_usage((self.find_alignments_dc, (s1, s2, score_matrix, local), ),
                                               retval=True, max_iterations=1)
        toc_dc = time.perf_counter()
        tic_dp = time.perf_counter()
        mem_usage_dp, dp_align = memory_usage((self.find_alignments_dp, (s1, s2, score_matrix, local), ),
                                               retval=True, max_iterations=1)
        toc_dp = time.perf_counter()
        if print_bool:
            print(f">> length {length} -> total length {length_sum}")
            print(f"   dp memory: {int(max(mem_usage_dp))} time: {toc_dp-tic_dp:0.4f}")
            print(f"   dc memory: {int(max(mem_usage_dc))} time: {toc_dc-tic_dc:0.4f}")
        
        num_aligns = len(dp_align)
        return toc_dp-tic_dp, toc_dc-tic_dc, max(mem_usage_dp), max(mem_usage_dc), num_aligns

    def alignment_cases(self, length, case=1, score_case=0, local=False, seed=1):
        # String cases
        if case == 0:
            s1_rand = self.random_string(int(length/4), seed=seed)
            s1_side_l = int(length * (3/4)/2)
            s1 = "A"*s1_side_l + s1_rand + "A"*s1_side_l
            s2 = s1

        if case == 1:
            s1_half = int(length/2)
            s1 = "A"*s1_half + "GTC" + "A"*s1_half
            s2 = "GTC"
        
        if case == 2:
            s1 = self.random_string(length, seed=seed)
            s2 = self.random_string(length, seed=seed+2)

        if case == 3:
            s1_half = int(length/2)
            s1 = "A"*s1_half + "GTC" + "A"*s1_half
            s2 = s1
        
        if case == 4:
            s1_quarter = int(length/4)
            s1 = "A"*s1_quarter + "T"*s1_quarter + "G"*s1_quarter + "C"*s1_quarter
            s2=s1
        
        # aligning sequence with itself, only that indels and snps are added.
        if case in [5,6,7]:
            s1 = self.random_string(length)
            s2 = ""
            random.seed(seed)
            s2_mismatch_inds = random.sample(range(len(s1)), int(length*0.1))
            random.seed(seed+5)
            s2_insert_inds = random.sample(range(len(s1)), int(length*0.1))
            for i in range(len(s1)):
                if case != 6 and i in s2_mismatch_inds:
                    s = {'A', 'C', 'G', 'T'}
                    s2 += random.choice(list(s))
                elif case != 7 and i in s2_insert_inds:
                    s2 += s1[i]
                    s = {'A', 'C', 'G', 'T'}
                    s2 += random.choice(list(s))
                else:
                    s2 += s1[i]

        # score cases
        if score_case ==1:
            score_matrix = self.create_score_matrix(match=1, mismatch=-1, indel=-1)
        if score_case == 2:
            score_matrix = self.create_score_matrix(match=10, mismatch=-15, indel=-15)
        if score_case == 3:
            score_matrix = self.create_score_matrix(match=10, mismatch=-15, indel=-7)
        if score_case == 4:
            score_matrix = self.create_score_matrix(match=10, mismatch=-8, indel=-15)
        
        
        return s1, s2, score_matrix

    def averages(self, case=0, score_case=2, lengths=[100], local=False, custom_str=""):
        loc_str = ""
        if local:
            loc_str = " for local alignment"
        else:
            loc_str = " for global alignment"
        print(f"\n\nAverages testing{loc_str}. {custom_str}")
        print("Press `ENTER` to proceed OR enter any key to skip.")
        q = input()
        if q.lower() != "":
            return
        dp_times = []; dc_times = []; dp_mems = []; dc_mems = []; num_aligns = []
        for l in lengths:
            for seed in range(1,4):
                s1, s2, score_matrix = self.alignment_cases(l, case=case, score_case=score_case, seed=seed)
                dp_t, dc_t, dp_m, dc_m , n_al= self.time_mem_sanity(s1, s2, score_matrix, length=l, local=local, print_bool=False)
                dp_times.append(dp_t)
                dc_times.append(dc_t)
                dp_mems.append(dp_m)
                dc_mems.append(dc_m)
                num_aligns.append(n_al)
            print(f">> 3*2 alignments for DP, DC: total length {len(s1)+len(s2)}")
            print(f"   dp time mean: {mean(dp_times):0.4f}")
            print(f"   dc time mean: {mean(dc_times):0.4f}")
            print(f"   dp mem mean: {mean(dp_mems):0.4f}")
            print(f"   dc mem mean: {mean(dc_mems):0.4f}")
            print(f"   Returned an average {int(mean(num_aligns))} alignment(s).")

    def test_a(self):
        custom_str = "Fully random string test."
        lengths = [10, 50, 100]
        lengths2 = [50, 100, 200]
        self.averages(case=2, score_case=4, lengths=lengths, custom_str=custom_str)
        self.averages(case=2, score_case=4, lengths=lengths2, local=True, custom_str=custom_str)

    def test_b(self):
        custom_str = "s1 = random string, s2 = s1 with added SNPs and indels."
        lengths = [10, 200, 300]
        self.averages(case=5, score_case=2, lengths=lengths, custom_str=custom_str)
        self.averages(case=5, score_case=2, lengths=lengths, local=True, custom_str=custom_str)

    def test_b_i(self):
        custom_str = "s1 = random string, s2 = s1 with added indels."
        lengths = [10, 200, 300]
        self.averages(case=6, score_case=3, lengths=lengths, custom_str=custom_str)
        self.averages(case=6, score_case=3, lengths=lengths, local=True, custom_str=custom_str)

    def test_b_ii(self):
        custom_str = "s1 = random string, s2 = s1 with added SNPs."
        lengths = [10, 100, 150]
        self.averages(case=7, score_case=4, lengths=lengths, custom_str=custom_str)
        self.averages(case=7, score_case=3, lengths=lengths, custom_str=custom_str)
        self.averages(case=7, score_case=4, lengths=lengths, local=True, custom_str=custom_str)
        self.averages(case=7, score_case=2, lengths=lengths, local=True, custom_str=custom_str)

    def test_c(self):
        custom_str = "Very long-short strings."
        lengths = [50000, 100000, 250000]
        self.averages(case=1, score_case=2, lengths=lengths, custom_str=custom_str)
        self.averages(case=1, score_case=2, lengths=lengths, local=True, custom_str=custom_str)

    def test_d(self):
        custom_str = "Long semi-random. THIS WILL TAKE A WHILE."
        lengths = [1000, 2000, 3000, 4000]
        self.averages(case=0, score_case=4, lengths=lengths, custom_str=custom_str)
        self.averages(case=0, score_case=4, lengths=lengths, local=True, custom_str=custom_str)

if __name__ == "__main__":
    unittest.main(verbosity=2)