def readGenome(filename):
    with open(filename) as f:
        f.readline()
        genome = "".join([line.strip() for line in f.readlines()])
        return genome
