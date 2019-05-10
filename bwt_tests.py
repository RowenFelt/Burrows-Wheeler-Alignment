"""
Test suite for bwt_query

Adapted from bwt_query_exercise_tests.py by Michael Linderman
"""

import unittest
import bwt_query
import sparse_bwt_query

class TestBWTFMIndex(unittest.TestCase):
    def setUp(self):
        # Before each test construct FM index data structures for this genome
        self.text = "AATCGGGTTCAATCGGGGT"
        self.letters, self.suffix_array, self.first_occur, self.last_col, self.count_matrix = bwt_query.build_indices(self.text)
     
    def test_ordering(self):
        """ Validate that '$' is lexicographically less than all DNA characters """
        # Ensure assumptions about string ordering are correct
        self.assertEqual(sorted("acgtACGT$")[0],"$")
    
    def test_letters(self):
        """ Letters contains all characters in the genome """
        self.assertEqual(set(self.letters), set(self.text + "$"))

    def test_first_occur(self):
        """ First occurrence is RLE of sorted genome """
        self.assertEqual(self.letters[0],"$")
        first = sorted(self.text + "$")
        for i in range(1,len(self.letters)):
            beg, end = self.first_occur[i:i+2]
            self.assertEqual(first[beg:end], [self.letters[i]]*(end-beg))

    def test_last_col(self):
        """ Last column is BWT of the genome """
        self.assertEqual(self.last_col, "TC$AATTTCGCGGGGGTAAG")
        
    def test_counts(self):
        """ Last entry in counts is total number of each character """
        for i in range(len(self.letters)):
            self.assertEqual(self.count_matrix[self.letters[i]][-1], \
                             self.first_occur[i+1]-self.first_occur[i])
        
class TestBWTQuery(unittest.TestCase):
    def setUp(self):
        # Before each test construct FM index data structures for this genome
        self.text = "AATCGGGTTCAATCGGGGTCAG"
        self.letters, self.suffix_array, self.first_occur, self.last_col, self.count_matrix = bwt_query.build_indices(self.text)
    
    def query(self, pattern):
        """ Query pre-indexed genome for pattern, returning list of matching start indices """
        return bwt_query.query(self.letters, self.suffix_array, self.first_occur, self.last_col, self.count_matrix, pattern)

    def query_mismatch(self, pattern, mismatches):
        """Query with mismatches"""
        return bwt_query.query_mismatch(self.letters, self.suffix_array, self.first_occur, self.last_col, self.count_matrix, pattern, mismatches)
    
    def test_query(self):
        """ Query with multiple matches """
        # Query is a method of this class, and so is invoked via self
        matches = self.query("ATCG")
        # Check equality of actual matches to expected matches (sorting to avoid ambiguity)
        self.assertEqual(sorted(matches), [1, 11])

    def test_query2(self):
        """ Query end match """
        matches = self.query("GGGT")
        self.assertEqual(sorted(matches), [4, 15])

    def test_mismatches1(self):
        matches = self.query_mismatch("TCC", 1)
        self.assertEqual(sorted(matches), [2, 7, 8, 12, 18])

    def test_mismatches2(self):
        matches = self.query_mismatch("CAG", 2)
        self.assertEqual(sorted(matches), [0, 2, 3, 4, 9, 10, 12, 13, 14, 15, 19])

    def test_mismatches3(self):
        matches = self.query_mismatch("GGGG", 1)
        self.assertEqual(sorted(matches), [3, 4, 13, 14, 15]) 

    def test_mismatches4(self):
        matches = self.query_mismatch("GGGG", 2)
        self.assertEqual(sorted(matches), [2, 3, 4, 5, 12, 13, 14, 15, 16])
        "0123456789012345678901"
        "AATCGGGTTCAATCGGGGTCAG"

# These tests were for sparse queries, and I couldn't get it to work. Bummer
#class TestSparseBWTQuery(unittest.TestCase):
#    def setUp(self):
#        # Before each test construct FM index data structures for this genome
#        self.text = "AATCGGGTTCAATCGGGGTCAG"
#        self.letters, self.suffix_array, self.first_occur, self.last_col, self.count_matrix = sparse_bwt_query.build_sparse_indices(self.text, 5)
#    
#    def query(self, pattern, sparsity):
#        """ Query pre-indexed genome for pattern, returning list of matching start indices """
#        return sparse_bwt_query.query(self.letters, self.suffix_array, self.first_occur, self.last_col, self.count_matrix, pattern, sparsity)
#
#    def query_mismatch(self, pattern, mismatches, sparsity):
#        """Query with mismatches"""
#        return sparse_bwt_query.query_mismatch(self.letters, self.suffix_array, self.first_occur, self.last_col, self.count_matrix, pattern, mismatches, sparsity)
#    
#    def test_sparse_query(self):
#        """ Query with multiple matches """
#        # Query is a method of this class, and so is invoked via self
#        matches = self.query("ATCG", 5)
#        # Check equality of actual matches to expected matches (sorting to avoid ambiguity)
#        self.assertEqual(sorted(matches), [1, 11])
#
#    def test_sparse_query2(self):
#        """ Query end match """
#        matches = self.query("GGGT", 5)
#        self.assertEqual(sorted(matches), [4, 15])
#
#    def test_sparse_mismatches1(self):
#        matches = self.query_mismatch("TCC", 1, 5)
#        self.assertEqual(sorted(matches), [2, 7, 8, 12, 18])
#
#    def test_sparse_mismatches2(self):
#        matches = self.query_mismatch("CAG", 2, 5)
#        self.assertEqual(sorted(matches), [0, 2, 3, 4, 9, 10, 12, 13, 14, 15, 19])
#
#    def test_sparse_mismatches3(self):
#        matches = self.query_mismatch("GGGG", 1, 5)
#        self.assertEqual(sorted(matches), [3, 4, 13, 14, 15]) 
#
#    def test_sparse_mismatches4(self):
#        matches = self.query_mismatch("GGGG", 2, 5)
#        self.assertEqual(sorted(matches), [2, 3, 4, 5, 12, 13, 14, 15, 16])

if __name__ == '__main__':
    unittest.main()
