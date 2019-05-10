"""
based on 'Rosalind.info problem 9N in-class exercise' written by Michael Linderman
"""

import sys
import bwt_sort
import bwt_read

def read9N(filename):
    """ Read 9N input """
    with open(filename, "r") as file:
        text = file.readline().strip()
        patterns = [line.strip() for line in file]
    return text, patterns

def build_first_occur(first_column):
    """ Build first occurrence array from first_column
    
    Args:
        first_column: First column as string
        
    Return:
        Tuple of (letters as string, list of first occurrences)
    """
    letters = "$"
    first_occur = [0]
    for i in range(1, len(first_column)):
        letter = first_column[i]
        if letter != letters[-1]:
            letters += letter
            first_occur.append(i)
    first_occur.append(len(first_column)) 
    
    return (letters, first_occur)

def build_count_matrix(letters, last_column):
    """ Build count matrix from last_column
    
    Args:
        letters: String of letters in index
        last_column: Last column as string
        
    Return:
        Dictionary of count array indexed by symbol
    """ 
    count_matrix = {}
    for letter in letters:
        count_matrix[letter] = [0]
    for current_letter in last_column:
        for letter in letters:
            letter_list = count_matrix[letter]
            if letter == current_letter:
                letter_list.append(letter_list[-1] + 1)
            else:
                letter_list.append(letter_list[-1])    
    return count_matrix

def build_sparse_count_matrix(letters, last_column, sparsity):
    """ Build count matrix from last_column
    
    Args:
        letters: String of letters in index
        last_column: Last column as string
        sparsity: factor by which to reduce count_matrix
        
    Return:
        Dictionary of count array indexed by symbol with sparsity
    """ 
    count_matrix = {}
    count_letters = {}
    for letter in letters:
        count_matrix[letter] = []
        count_letters[letter] = 0
    i = 0
    for current_letter in last_column:
        if i % sparsity == 0:
            for letter in letters:
                count_matrix[letter].append(count_letters[letter])
        i+=1
        count_letters[current_letter] = count_letters[current_letter] + 1 #sum of letters
    return count_matrix

def build_sparse_indices(text, sparsity):
    """Build BWT FM-index data structures
    
    Args:
        text: String to index
    
    Returns:
        Tuple of
        (letters, suffix array, first occurrence, last col, count matrix)
    """
    assert(sorted(text + "$")[0] == "$")
    
    # Append "$"
    text = text + "$"
    G = len(text)
    
    # Create all cyclic rotations with starting position
    cycles = [(index, (G-1+index) % G) for index in range(G)]
    
    # sort cycles lexographically with pointer based quicksort
    bwt_sort.quickSort(cycles, 0, G-1, text, G)
      
    suffix_array = [cycle[0] for cycle in cycles] 
    last_col = "".join([text[cycle[1]] for cycle in cycles])   
    
    # Determine 1st occurrence
    first_col = "".join([text[cycle[0]] for cycle in cycles])
    letters, first_occur = build_first_occur(first_col)

    count_matrix = build_sparse_count_matrix(letters, last_col, sparsity)
        
    return (letters, suffix_array, first_occur, last_col, count_matrix)


def query(letters, suffix_array, first_occur, last_col, count_matrix, pattern, sparsity):
    """ A wrapper for querying without mismatches """
    return query_bwt(letters, suffix_array, first_occur, last_col, count_matrix, pattern, 0, None, None, sparsity)

def query_mismatch(letters, suffix_array, first_occur, last_col, count_matrix, pattern, mismatches, sparsity):
    """ A wrapper for querying with mismatches """
    matches = set()
    base = query_bwt(letters, suffix_array, first_occur, last_col, count_matrix, pattern, mismatches, None, None, sparsity) 
    for match in base:
        matches.add(match)
    for i in range(1,mismatches+1):
        new_pattern = pattern[:-i]
        temp = query_bwt(letters, suffix_array, first_occur, last_col, count_matrix, new_pattern, mismatches-i, None, None, sparsity)
        for match in temp:
            matches.add(match)
    return matches
    

def query_bwt(letters, suffix_array, first_occur, last_col, count_matrix, pattern, mismatches, top, bot, sparsity):   
    """ Query BWT FM-index
    
    Args:
        letters: String of letters in index
        suffix_array: SuffixArray as list
        first_occur: List of first occurrences
        last_col: Last column of BW matrix
        count_matrix: Counts as a dictionary
        pattern: Query as a string
        mismatches: number of allowed mismatches
        top: top of current view of BWT in recursive call, None if first call
        bot: bottom of current view of BWT in recursive call, None if first call
        sparsity: sparsity factor of checkpoint matrix
    
    Returns:
        List of match positions given a number of mismatches, last letter cannot be a mismatch
    """

    # Find initial top and bottom pointers
    if top == None or bot == None:
        which_letter = letters.find(pattern[-1])
        top = first_occur[which_letter]
        bot = first_occur[which_letter+1]    
    
    matches = set()
    i = len(pattern)-2
    while i >= 0 and top != bot:
        if mismatches > 0:
            for letter in letters:
                if letter == "$": #if it's the $ skip it
                    continue
                new_top, new_bot = calculate_sparse_count(letters, letter, i, first_occur, count_matrix, top, bot, last_col, sparsity)
                if new_top != None and new_bot != None:
                    new_pattern = pattern[:i] + letter          #create new pattern with letter 
                    if letter == pattern[i]:
                        new_matches = query_bwt(letters, suffix_array, first_occur, last_col, count_matrix, new_pattern, mismatches, new_top, new_bot, sparsity)
                    else:
                        new_matches = query_bwt(letters, suffix_array, first_occur, last_col, count_matrix, new_pattern, mismatches-1, new_top, new_bot, sparsity)
                    for match in new_matches:
                        matches.add(match)
            return matches
        else: #no more mismatches
            top, bot = calculate_sparse_count(letters, pattern[i], i, first_occur, count_matrix, top, bot, last_col, sparsity)
            i -= 1

    for suffix in suffix_array[top:bot]:        #add direct matches
        matches.add(suffix)
    return list(matches)

def calculate_sparse_count(letters, letter, index, first_occur, count_matrix, top, bot, last_col, sparsity):
    """ Calculates the new top and bottom pointers using the sparse count matrix 
    
    Args: 
        letters: String of letters in index
        suffix_array: SuffixArray as list
        first_occur: List of first occurrences
        last_col: Last column of BW matrix as string
        count_matrix: Counts as a dictionary
        top: top of current view of BWT in recursive call, None if first call
        bot: bottom of current view of BWT in recursive call, None if first call
        sparsity: sparsity factor of checkpoint matrix
    
    Returns:
        new top and bottom pointers
        returns None, None if the given letter does not appear in the range
    """
    first = first_occur[letters.find(letter)] # index of first occurrence in first_col
    # Generate "next" top and bottom pointers using count array
    count_array = count_matrix[letter] #count array for a given symbol
    top_const = count_array[top//sparsity] #find most recently populated row
                                        #analagous to count_array[bot] in non-sparse count matrix
    for i in range(1 + (top//sparsity) * sparsity, top+1):
        if last_col[i-1] == letter:
            top_const += 1   
    bot_const = count_array[bot//sparsity] #analagous to count_array[bot] in non-sparse count matrix
    for i in range(1 + (bot//sparsity) * sparsity, bot+1):
        if last_col[i-1] == letter:
            bot_const += 1   
    if bot_const - top_const == 0: #if the letter does not appear in this range argv
        # might also need something like x and 
        return None, None
    top = first + top_const
    bot = first + bot_const
    return top, bot

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(f"Usage: python3 {__file__} INPUT_FILE")
        sys.exit(1)

    text = bwt_read.readGenome(sys.argv[1])
    letters, suffix_array, first_occur, last_col, count_matrix = build_indices(text)
    matches = []
    for pattern in patterns:
        matches += query(letters, suffix_array, first_occur, last_col, count_matrix, pattern)
    print(" ".join(map(str, sorted(matches))))
        
