"""
    An experimental setup.
    Takes two inputs:
        argv[1]:    a file containing the reference genome
        argv[2]:    a file containing the experimental reads
        the following experimental results will be printed to the terminal for each run:
            how many reads resulted in a match
            how many reads matched the correct position
            how many aver
"""

import bwt_query
import bwt_read
import sys
import time


if __name__ == "__main__":
    if len(sys.argv) < 4: 
        print("Usage: experiment.py <reference-genome> <experimental-reads> <query_results_file(output)>")
        sys.exit(1) 

    text = bwt_read.readGenome(sys.argv[1]) 
    queries = bwt_read.readSAM(sys.argv[2])
    letters, suffix_array, first_occur, last_col, count_matrix = bwt_query.build_indices(text) 
    for i in range(0,5):
        total_matches = 0
        matches_to_genome = 0
        correct_queries = 0
        incorrect_matches = 0
        start = time.time()
        for query in queries:
            matches = bwt_query.query_mismatch(letters, suffix_array, first_occur, last_col, count_matrix, query[3], i)
            if matches:
                matches_to_genome += 1
            for match in matches:
                if match+1 == int(query[1]):
                    correct_queries += 1
                else:
                    incorrect_matches += 1
                total_matches += 1
        end = time.time()
        average_matches_per_read = total_matches / len(queries)
        print("For " + str(i)+" mismatches:")
        print("\t" + str(matches_to_genome) + " reads matched to the genome")
        print("\t" + str(correct_queries) + " matches mapped to the correct position")
        print("\t" + str(incorrect_matches) + " incorrect mappings")
        print("\t" + str(average_matches_per_read) + " average matches per read")
        print("\tExecuted in " + str(end - start) + " seconds")

