def readGenome(filename):
    with open(filename) as f:
        f.readline()
        genome = "".join([line.strip() for line in f.readlines()])
        return genome

def readSAM(filename):
    """ Returns a list of all the queries formatted as
        [id, index, cigar, query string]
    """
    with open(filename) as f:
        f.readline()
        f.readline()
        f.readline()
        queries = []
        for line in f.readlines():
            line = line.split()
            queries.append([line[2], line[3], line[5], line[9]])
        return queries
