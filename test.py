import unittest
import os, sys

def sample(input_filename, m, n, fn, fp):
    command = f'build/src-unigen/lumberjack {input_filename} -c {m} -m {n} -n {fn} -p {fp} | grep -Pzo \'.*solutions sampled(.*\\n)*\''
    os.system(command)

class TestDolloSAT(unittest.TestCase):

    # A single forbidden matrix, no clustering allowed, no false positives or false negatives allowed,
    # all losses allowed.
    #
    # 3 solutions are:
    # 1 2   1 0   1 2
    # 0 1   2 1   2 1
    # 1 1 , 1 1,  1 1
    def test_simple_allow_all_losses(self):
        print('=================== test_simple_allow_all_losses ===================')
        sample('test_inputs/simple_forbidden.txt', 3, 2, 0, 0)
        print('Expected # solutions: 3')

    # No way to cluster this matrix to a 4x4 using only one false negative.
    #
    # No solutions
    def test_no_solutions_no_clustering(self):
        print('=================== test_no_solutions_no_clustering ===================')
        print('Expected # solutions: 0')
        # false neg rate = 1 / # zeroes = 1 / 12 
        sample('test_inputs/no_clustering.txt', 4, 4, 0.084, 0)
    
    # No false positives or false negatives allowed.
    #
    # 6 possible solutions for given matrix:
    # 1 0 0  1 2 0  1 2 2  1 2 2  1 0 0  1 2 0  
    # 1 1 0  1 1 0  1 1 0  1 1 2  1 1 2  1 1 2  
    # 1 1 1, 1 1 1, 1 1 1, 1 1 1, 1 1 1, 1 1 1
    def test_harder_no_error(self):
        print('=================== test_harder_no_error ===================')
        print('Expected # solutions: 6')
        sample('test_inputs/test_harder.txt', 3, 3, 0, 0)
    
    # One false negative
    # 
    # False negative at B[0][1] -> clusters to a 2 x 2 matrix, no solutions
    # False negative at B[0][2] -> 3 solutions
    # False negative at B[1][2] -> clusters to a 2x2 matrix, no solutions
    #
    # Total = 6 + 3 = 9 solutions
    def test_harder_one_fn(self):
        print('=================== test_harder_one_fn ===================')
        print('Expected # solutions: 3')
        # false neg rate = 1 / # zeroes = 1 / 3
        sample('test_inputs/test_harder.txt', 3, 3, 0.334, 0)
    
    # A matrix that clusters to a forbidden matrix, no false positives or
    # false negatives allowed, all mutation loss is allowed.
    # 
    # 3 solutions (should be same as test_simple_allow_all_losses)
    def test_cell_cluster_to_forbidden_allow_losses(self):
        print('=================== test_cell_cluster_to_forbidden_allow_losses ===================')
        print('Expected # solutions: 3')
        sample('test_inputs/cluster_cells.txt', 3, 2, 0, 0)
    
    # A matrix that clusters to a forbidden matrix, no false positives or
    # false negatives allowed, all mutation loss is allowed.
    # 
    # 3 solutions (should be same as test_simple_allow_all_losses)
    def test_mutation_cluster_to_forbidden_allow_losses(self):
        print('=================== test_mutation_cluster_to_forbidden_allow_losses ===================')
        print('Expected # solutions: 3')
        sample('test_inputs/cluster_mutations.txt', 3, 2, 0, 0)
    
    # Matrix should cluster to this matrix, with no false negatives or positives allowed:
    # 1 0 0
    # 1 1 0
    # 1 1 1
    #
    # 6 possible solutions for this clustered matrix:
    # 1 0 0  1 2 0  1 2 2  1 2 2  1 0 0  1 2 0  
    # 1 1 0  1 1 0  1 1 0  1 1 2  1 1 2  1 1 2  
    # 1 1 1, 1 1 1, 1 1 1, 1 1 1, 1 1 1, 1 1 1
    def test_harder_cluster(self):
        print('=================== test_harder_cluster ===================')
        print('Expected # solutions: 6')
        sample('test_inputs/test_harder_clustering.txt', 3, 3, 0, 0)
    
    # Input matrix of all 1s, 2 false positives allowed, 0 false negatives, all allowed losses
    #
    # All row permutations of the following 5 solutions are valid:
    # 1 2   1 0   1 0   0 1   0 1   
    # 2 1   2 1   1 2   2 1   1 2  
    # 1 1 , 1 1 , 1 1 , 1 1 , 1 1
    # 5*3! = 30 solutions
    def test_more_fp(self):
        print('=================== test_more_fp ===================')
        print('Expected # solutions: 30')
        # false pos rate = 2 / # ones = 2 / 6
        sample('test_inputs/ones_3x2.txt', 3, 2, 0, 0.334)
    
    # Input matrix is 4x4, clustered to 2x2 matrix, 2 false positives, 0 false negatives, all losses
    #
    # 2 solutions:
    # 1 0  1 2
    # 1 1, 1 1
    def test_small_cluster(self):
        print('=================== test_small_cluster ===================')
        print('Expected # solutions: 2')
        # false pos rate = 1 / # ones = 1 / 11
        sample('test_inputs/cluster_small.txt', 2, 2, 0, 0.091)

if __name__ == '__main__':
    unittest.main()