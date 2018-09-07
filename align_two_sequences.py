#!/usr/bin/python

"""
Author: Elena Garcia
Course: Algorithms in Bioinformatics
"""

# for commandline options:
from optparse import OptionParser, OptionGroup

# Built-in exchange matrices.
identity = [
    [ 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,],
    [-1, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,],
    [-1,-1, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,],
    [-1,-1,-1, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,],
    [-1,-1,-1,-1, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,],
    [-1,-1,-1,-1,-1, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,],
    [-1,-1,-1,-1,-1,-1, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,],
    [-1,-1,-1,-1,-1,-1,-1, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,],
    [-1,-1,-1,-1,-1,-1,-1,-1, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,],
    [-1,-1,-1,-1,-1,-1,-1,-1,-1, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,],
    [-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,],
    [-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,],
    [-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,],
    [-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,],
    [-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,],
    [-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,],
    [-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,],
    [-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 1,-1,-1,-1,-1,-1,-1,-1,-1,],
    [-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 1,-1,-1,-1,-1,-1,-1,-1,],
    [-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 1,-1,-1,-1,-1,-1,-1,],
    [-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 1,-1,-1,-1,-1,-1,],
    [-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 1,-1,-1,-1,-1,],
    [-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 1,-1,-1,-1,],
    [-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 1,-1,-1,],
    [-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 1,-1,],
    [-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 1,],
]

pam250 = [
    [ 2, 0,-2, 0, 0,-4, 1,-1,-1, 0,-1,-2,-1, 0, 0, 1, 0,-2, 1, 1, 0, 0,-6, 0,-3, 0,],
    [ 0, 2,-4, 3, 2,-5, 0, 1,-2, 0, 1,-3,-2, 2, 0,-1, 1,-1, 0, 0, 0,-2,-5, 0,-3, 2,],
    [-2,-4,12,-5,-5,-4,-3,-3,-2, 0,-5,-6,-5,-4, 0,-3,-5,-4, 0,-2, 0,-2,-8, 0, 0,-5,],
    [ 0, 3,-5, 4, 3,-6, 1, 1,-2, 0, 0,-4,-3, 2, 0,-1, 2,-1, 0, 0, 0,-2,-7, 0,-4, 3,],
    [ 0, 2,-5, 3, 4,-5, 0, 1,-2, 0, 0,-3,-2, 1, 0,-1, 2,-1, 0, 0, 0,-2,-7, 0,-4, 3,],
    [-4,-5,-4,-6,-5, 9,-5,-2, 1, 0,-5, 2, 0,-4, 0,-5,-5,-4,-3,-3, 0,-1, 0, 0, 7,-5,],
    [ 1, 0,-3, 1, 0,-5, 5,-2,-3, 0,-2,-4,-3, 0, 0,-1,-1,-3, 1, 0, 0,-1,-7, 0,-5,-1,],
    [-1, 1,-3, 1, 1,-2,-2, 6,-2, 0, 0,-2,-2, 2, 0, 0, 3, 2,-1,-1, 0,-2,-3, 0, 0, 2,],
    [-1,-2,-2,-2,-2, 1,-3,-2, 5, 0,-2, 2, 2,-2, 0,-2,-2,-2,-1, 0, 0, 4,-5, 0,-1,-2,],
    [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,],
    [-1, 1,-5, 0, 0,-5,-2, 0,-2, 0, 5,-3, 0, 1, 0,-1, 1, 3, 0, 0, 0,-2,-3, 0,-4, 0,],
    [-2,-3,-6,-4,-3, 2,-4,-2, 2, 0,-3, 6, 4,-3, 0,-3,-2,-3,-3,-2, 0, 2,-2, 0,-1,-3,],
    [-1,-2,-5,-3,-2, 0,-3,-2, 2, 0, 0, 4, 6,-2, 0,-2,-1, 0,-2,-1, 0, 2,-4, 0,-2,-2,],
    [ 0, 2,-4, 2, 1,-4, 0, 2,-2, 0, 1,-3,-2, 2, 0,-1, 1, 0, 1, 0, 0,-2,-4, 0,-2, 1,],
    [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,],
    [ 1,-1,-3,-1,-1,-5,-1, 0,-2, 0,-1,-3,-2,-1, 0, 6, 0, 0, 1, 0, 0,-1,-6, 0,-5, 0,],
    [ 0, 1,-5, 2, 2,-5,-1, 3,-2, 0, 1,-2,-1, 1, 0, 0, 4, 1,-1,-1, 0,-2,-5, 0,-4, 3,],
    [-2,-1,-4,-1,-1,-4,-3, 2,-2, 0, 3,-3, 0, 0, 0, 0, 1, 6, 0,-1, 0,-2, 2, 0,-4, 0,],
    [ 1, 0, 0, 0, 0,-3, 1,-1,-1, 0, 0,-3,-2, 1, 0, 1,-1, 0, 2, 1, 0,-1,-2, 0,-3, 0,],
    [ 1, 0,-2, 0, 0,-3, 0,-1, 0, 0, 0,-2,-1, 0, 0, 0,-1,-1, 1, 3, 0, 0,-5, 0,-3,-1,],
    [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,],
    [ 0,-2,-2,-2,-2,-1,-1,-2, 4, 0,-2, 2, 2,-2, 0,-1,-2,-2,-1, 0, 0, 4,-6, 0,-2,-2,],
    [-6,-5,-8,-7,-7, 0,-7,-3,-5, 0,-3,-2,-4,-4, 0,-6,-5, 2,-2,-5, 0,-6,17, 0, 0,-6,],
    [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,],
    [-3,-3, 0,-4,-4, 7,-5, 0,-1, 0,-4,-1,-2,-2, 0,-5,-4,-4,-3,-3, 0,-2, 0, 0,10,-4,],
    [ 0, 2,-5, 3, 3,-5,-1, 2,-2, 0, 0,-3,-2, 1, 0, 0, 3, 0, 0,-1, 0,-2,-6, 0,-4, 3,],
]

blosum62 = [
    [ 4,-2, 0,-2,-1,-2, 0,-2,-1, 0,-1,-1,-1,-2, 0,-1,-1,-1, 1, 0, 0, 0,-3, 0,-2,-1,],
    [-2, 4,-3, 4, 1,-3,-1, 0,-3, 0, 0,-4,-3, 3, 0,-2, 0,-1, 0,-1, 0,-3,-4,-1,-3, 1,],
    [ 0,-3, 9,-3,-4,-2,-3,-3,-1, 0,-3,-1,-1,-3, 0,-3,-3,-3,-1,-1, 0,-1,-2,-2,-2,-3,],
    [-2, 4,-3, 6, 2,-3,-1,-1,-3, 0,-1,-4,-3, 1, 0,-1, 0,-2, 0,-1, 0,-3,-4,-1,-3, 1,],
    [-1, 1,-4, 2, 5,-3,-2, 0,-3, 0, 1,-3,-2, 0, 0,-1, 2, 0, 0,-1, 0,-2,-3,-1,-2, 4,],
    [-2,-3,-2,-3,-3, 6,-3,-1, 0, 0,-3, 0, 0,-3, 0,-4,-3,-3,-2,-2, 0,-1, 1,-1, 3,-3,],
    [ 0,-1,-3,-1,-2,-3, 6,-2,-4, 0,-2,-4,-3, 0, 0,-2,-2,-2, 0,-2, 0,-3,-2,-1,-3,-2,],
    [-2, 0,-3,-1, 0,-1,-2, 8,-3, 0,-1,-3,-2, 1, 0,-2, 0, 0,-1,-2, 0,-3,-2,-1, 2, 0,],
    [-1,-3,-1,-3,-3, 0,-4,-3, 4, 0,-3, 2, 1,-3, 0,-3,-3,-3,-2,-1, 0, 3,-3,-1,-1,-3,],
    [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,],
    [-1, 0,-3,-1, 1,-3,-2,-1,-3, 0, 5,-2,-1, 0, 0,-1, 1, 2, 0,-1, 0,-2,-3,-1,-2, 1,],
    [-1,-4,-1,-4,-3, 0,-4,-3, 2, 0,-2, 4, 2,-3, 0,-3,-2,-2,-2,-1, 0, 1,-2,-1,-1,-3,],
    [-1,-3,-1,-3,-2, 0,-3,-2, 1, 0,-1, 2, 5,-2, 0,-2, 0,-1,-1,-1, 0, 1,-1,-1,-1,-1,],
    [-2, 3,-3, 1, 0,-3, 0, 1,-3, 0, 0,-3,-2, 6, 0,-2, 0, 0, 1, 0, 0,-3,-4,-1,-2, 0,],
    [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,],
    [-1,-2,-3,-1,-1,-4,-2,-2,-3, 0,-1,-3,-2,-2, 0, 7,-1,-2,-1,-1, 0,-2,-4,-2,-3,-1,],
    [-1, 0,-3, 0, 2,-3,-2, 0,-3, 0, 1,-2, 0, 0, 0,-1, 5, 1, 0,-1, 0,-2,-2,-1,-1, 3,],
    [-1,-1,-3,-2, 0,-3,-2, 0,-3, 0, 2,-2,-1, 0, 0,-2, 1, 5,-1,-1, 0,-3,-3,-1,-2, 0,],
    [ 1, 0,-1, 0, 0,-2, 0,-1,-2, 0, 0,-2,-1, 1, 0,-1, 0,-1, 4, 1, 0,-2,-3, 0,-2, 0,],
    [ 0,-1,-1,-1,-1,-2,-2,-2,-1, 0,-1,-1,-1, 0, 0,-1,-1,-1, 1, 5, 0, 0,-2, 0,-2,-1,],
    [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,],
    [ 0,-3,-1,-3,-2,-1,-3,-3, 3, 0,-2, 1, 1,-3, 0,-2,-2,-3,-2, 0, 0, 4,-3,-1,-1,-2,],
    [-3,-4,-2,-4,-3, 1,-2,-2,-3, 0,-3,-2,-1,-4, 0,-4,-2,-3,-3,-2, 0,-3,11,-2, 2,-3,],
    [ 0,-1,-2,-1,-1,-1,-1,-1,-1, 0,-1,-1,-1,-1, 0,-2,-1,-1, 0, 0, 0,-1,-2,-1,-1,-1,],
    [-2,-3,-2,-3,-2, 3,-3, 2,-1, 0,-2,-1,-1,-2, 0,-3,-1,-2,-2,-2, 0,-1, 2,-1, 7,-2,],
    [-1, 1,-3, 1, 4,-3,-2, 0,-3, 0, 1,-3,-1, 0, 0,-1, 3, 0, 0,-1, 0,-2,-3,-1,-2, 4,],
]

def parse_commandline():
    usage = "%prog <fasta> [options]"
    version = "1.0"
    description = \
        "%prog aligns two sequences."
    epilog = \
        "Elena Garcia 2016"
    parser = OptionParser(usage=usage, description=description,
                          version="%prog "+version, epilog=epilog)

    # sequence/alignment options:
    parser.add_option("-v", dest="print_scorematrix", metavar="<file>",
                     help="input alignment file (fasta)")

    parser.add_option("-f", "--fasta",  dest="fasta", metavar="<file>",
                     help="input alignment file (fasta)")
    parser.set_defaults(fasta=None)

    parser.add_option("-e", "",  dest="exchange_matrix",
                     help="Exchange matrix: pam250, blosum62 or identity (%default%)")
    parser.set_defaults(exchange_matrix="pam250")

    parser.add_option("-l", "",  dest="align_local",  action="store_true",
                     help="align local")
    parser.set_defaults(align_local=False)

    parser.add_option("-g", "",  dest="align_global", action="store_true",
                     help="align global")
    parser.set_defaults(align_global=False)

    parser.add_option("-s", "",  dest="align_semiglobal", action="store_true",
                     help="align semi-global")
    parser.set_defaults(align_semiglobal=False)

    parser.add_option("-p", "",  dest="gap_penalty", type="int",
                     help="Gap penalty (%default%)")
    parser.set_defaults(gap_penalty=2)

    # get the options:
    (options, args) = parser.parse_args()

    if options.print_scorematrix:
        options.fasta = options.print_scorematrix

    if not options.fasta:
        # check if we have an option left (to be used as input filename):
        if args:
            options.fasta = args.pop()
        else:
            print "Need at least an input file (fasta)"
            print ""
            parser.print_help()
            print ""
            print "ERROR: no input file given"
            exit(-1)

    # check alignment type:
    align_options = [options.align_local, options.align_global, options.align_semiglobal]
    # check if at least one alignment option was true, else choose global
    if align_options.count(True)==0:
        print "No alignment type given, using Global"
        options.align_global=True
    # check if not more than one alignment option was true, else error and exit
    if align_options.count(True)>1:
        print "ERROR: multiple alignment types chosen"
        exit(-1)

    # check for any leftover command line arguments:
    if len(args):
        warning("ignoring additional arguments "+str(args))

    # clean up (recommended):
    del(parser)
    return options

class Sequence:
    """Stores a sequence object"""

    def __init__(self, Label="", Sequence="" ):
        """Initialize a new Sequence object

        Label -- identifier of sequence (text)
        Sequence -- sequence string in single-letter alphabet
        """
        self.Label       = Label
        self.Sequence    = Sequence

    # this makes that you can do 'print sequence' and get nice output:
    def __str__(self):
        """Return string representation of a Sequence object"""
        # newline-delimited values of all the attributes
        return ">%s\n%s" % (self.Label, self.Sequence)


def readSequences(lines):
    """Return Sequences object

    lines -- list of lines or any object that behaves like it

    This routine parses a fasta file and returns a list of Sequence objects
    containing the sequences with label and sequence data set
    """
    seqs = []
    label = None
    seq_lines = []
    for line in lines:
        line = line.strip()      # strip off white space
        if not line:             # skip empty lines
            continue
        if line.startswith(';'): # ignore comment lines
            continue
        # check for start of next sequence:
        if line.startswith('>'): # label line
            # first, store the previous sequence if we had one:
            if seq_lines:
                seqs.append(Sequence(label, ''.join(seq_lines)))
                seq_lines = []
            # get the label (name) for the next sequence
            label = line[1:].strip()
        else:
            # collect all lines with sequence information for this sequence:
            seq_lines.append(line)
    # take care of the last sequence in the file
    seqs.append(Sequence(label, ''.join(seq_lines)))
    return seqs

def do_global_alignment(sequences, matrix, penalty, alignment, printing):
    """ do pairwise GLOBAL alignment using DP
        At the end, it prints the alignment
    """
    #Gives easy name to sequences, and makes sure they're capital letters
    seq1 = sequences[0].Sequence
    seq2 = sequences[1].Sequence

    #Makes sure penalty is a negative value
    penalty = -(abs(penalty))

    # Creates the score matrix for the alignment
    score_matrix, traceback_matrix = create_score_matrix(seq1, seq2, matrix, penalty, alignment)

    # Recognizes the parse -v, and prints the score matrix only if -v is present
    if printing == True:
        print_matrix(seq1, seq2, score_matrix)

    # Aligns the sequences according to the score matrix
    i_max = len(seq1)-1
    j_max = len(seq2)-1
    aligned_sequences = traceback(seq1, seq2, score_matrix, traceback_matrix, i_max, j_max, alignment)
    print aligned_sequences

def do_semiglobal_alignment(sequences, matrix, penalty, alignment, printing):
    """ do pairwise SEMI-GLOBAL alignment using DP
        At the end, it prints the alignment
    """
    #Gives easy name to sequences, and makes sure they're capital letters
    seq1 = sequences[0].Sequence
    seq2 = sequences[1].Sequence

    #Makes sure penalty is a negative value
    penalty = -(abs(penalty))

    #Creates the score matrix for the alignment
    score_matrix, traceback_matrix = create_score_matrix(seq1, seq2, matrix, penalty, alignment)

    #Recognizes the parse -v, and prints the score matrix only if -v is present
    if printing == True:
        print_matrix(seq1, seq2, score_matrix)

    # Aligns the sequences according to the score matrix
    i_max, j_max, irrelevant = find_max_score(seq1, seq2, score_matrix, alignment)
    aligned_sequences = traceback(seq1, seq2, score_matrix, traceback_matrix, i_max, j_max, alignment)
    print aligned_sequences

def do_local_alignment(sequences, matrix, penalty, alignment, printing):
    """ do pairwise LOCAL alignment using DP
        At the end, it prints the alignment
    """
    #Gives easy name to sequences, and makes sure they're capital letters
    seq1 = sequences[0].Sequence
    seq2 = sequences[1].Sequence

    #Makes sure penalty is a negative value
    penalty = -(abs(penalty))

    #Creates the score matrix for the alignment
    score_matrix, traceback_matrix = create_score_matrix(seq1, seq2, matrix, penalty, alignment)

    #Recognizes the parse -v, and prints the score matrix only if -v is present
    if printing == True:
        print_matrix(seq1, seq2, score_matrix)

    # Finds the maximum score(s) across the matrix
    i_max, j_max, multiple_maxs = find_max_score(seq1, seq2, score_matrix, alignment)

    # Aligns the sequences according to the score matrix
    if multiple_maxs == True:
        aligned_sequences = ""
        n = 0
        for num in range(len(i_max)):
            maximum_i = i_max[num]
            maximum_j = j_max[num]
            n += 1
            one_aligned_sequence = traceback(seq1, seq2, score_matrix, traceback_matrix, maximum_i, maximum_j, alignment)
            aligned_sequences = aligned_sequences + "Alignment " + str(n) + ":\n" + one_aligned_sequence + "\n\n"
        aligned_sequences = aligned_sequences[:-2]
    else:
        aligned_sequences = traceback(seq1, seq2, score_matrix, traceback_matrix, i_max, j_max, alignment)

    # Prints the alignment(s)
    print aligned_sequences


def create_score_matrix(seq1, seq2, matrix, penalty, alignment):
    """Returns a score matrix
            Sequence 2 on top of the matrix, sequence 1 at the side

        It also returns a traceback matrix that will indicate the way-back
            in the traceback process
    """

    # Creates a matrix as a list of lists, that is filled with zeros
    columns = len(seq2)+1
    rows = len(seq1)+1
    score_matrix = [([0]*columns) for x in xrange(0,rows)]
    traceback_matrix = [([0]*columns) for x in xrange(0,rows)]

    # Computes entries in the score matrix and in the traceback matrix
    for i in range(rows):
        for j in range(columns):
            if alignment == "global":
                score_matrix[0][j] = j*penalty #First row filled with "0, -2, -4, ..."
                score_matrix[i][0] = i*penalty #First column filled with "0, -2, -4, ..."

            match = 0
            gap1 = 0
            gap2 = 0
            zero = 0

            #Score for match or mismatch
            aa1 = seq1[i-1]
            aa2 = seq2[j-1]
            if i>=1 and j>=1:
                match = score_matrix[i-1][j-1] + matrix[ord(aa1)-ord("A")][ord(aa2)-ord("A")]

            #Score for gap in seq1 (add -gap from left-cell)
            if j>=1:
                gap1 = score_matrix[i][j-1] + penalty

            #Score for gap in seq2 (add -gap from up-cell)
            if i>=1:
                gap2 = score_matrix[i-1][j] + penalty

            #Each cell gets the maximum score
            #The traceback matrix is filled
            if alignment in ["global", "semi"]:
                score_matrix[i][j] = max(match, gap1, gap2)
                if max(match, gap1, gap2) == match:
                    traceback_matrix[i][j] = "diag"
                if max(match, gap1, gap2) == gap1:
                    traceback_matrix[i][j] = "left"
                if max(match, gap1, gap2) == gap2:
                    traceback_matrix[i][j] = "up"

            if alignment == "local":
                score_matrix[i][j] = max(match, gap1, gap2, zero)
                if max(match, gap1, gap2, zero) == match:
                    traceback_matrix[i][j] = "diag"
                if max(match, gap1, gap2, zero) == gap1:
                    traceback_matrix[i][j] = "left"
                if max(match, gap1, gap2, zero) == gap2:
                    traceback_matrix[i][j] = "up"
                elif max(match, gap1, gap2, zero) == zero:
                    traceback_matrix[i][j] = "zero"

    return score_matrix, traceback_matrix

def find_max_score(seq1, seq2, score_matrix, alignment):
    """ In semi-global and local, finds the maximum score_matrix
        from where to start the traceback
    """

    rows = len(seq1)
    columns = len(seq2)
    i_max = len(seq1)
    j_max = len(seq2)

    multiple_maxs = False

    # Finds the highest score in the last column and last row
    if alignment == "semi":
        for i in range(0,rows):
            for j in range(0,columns):
                if score_matrix[i][columns] >= score_matrix[i_max][columns]:
                        i_max = i
                if score_matrix[rows][j] >= score_matrix[rows][j_max]:
                        j_max = j

        # Chooses the down-most and right-most possible maximum value
        if score_matrix[i_max][columns] > score_matrix[rows][j_max]:
            i = i_max-1
            j = columns-1
        elif score_matrix[i_max][columns] < score_matrix[rows][j_max]:
            j = j_max-1
            i = rows-1
        else:
            if i_max > j_max:
                i = i_max-1
                j = columns-1
            else:
                j = j_max-1
                i = rows-1

    # Finds the highest score across all the score matrix
    # First, it finds the right-most down-most highest score
    if alignment == "local":
        for i in range(0, rows):
            for j in range(0, columns):
                    if score_matrix[i][j_max] >= score_matrix[i_max][j_max]:
                        i_max = i
                    if score_matrix[i_max][j] >= score_matrix[i_max][j_max]:
                        j_max = j
        i = i_max-1
        j = j_max-1

    # Second, it finds other possible maximums and appends their position to the list
        i_list = [i]
        j_list = [j]
        if score_matrix[i_max][j_max] > 0:
            for item in range(0, i_max):
                for jtem in range(0, j_max):
                    if score_matrix[item][jtem] == score_matrix[i_max][j_max]:
                        multiple_maxs = True
                        i_list.append(item-1)
                        j_list.append(jtem-1)

    if multiple_maxs == True:
        return i_list, j_list, multiple_maxs
    else:
        return i, j, None

def traceback(seq1, seq2, score_matrix, traceback_matrix, i_max, j_max, alignment):
    """ Performs the traceback,
    returns aligned sequences """
    # Allocates empty strings for the aligned sequences
    seq1_traceback = ""
    seq2_traceback = ""
    matches_traceback = ""

    i = i_max+1
    j = j_max+1

    """ The traceback itselt: """
    # Writes the aligned sequences, taking into account the traceback matrix.
    # It adds gaps when necessary.
    # It displays the alignment with "|" when mathces occur.
    for n in range(0,(i+1)*(j+1)):
        if i>0 and j>0:
            if traceback_matrix[i][j] == "zero":
                break
            if traceback_matrix[i][j] == "diag":
                seq1_traceback = seq1[i-1] + seq1_traceback
                seq2_traceback = seq2[j-1] + seq2_traceback
                if seq1[i-1] == seq2[j-1]:
                    matches_traceback = "|" + matches_traceback
                else:
                    matches_traceback = " " + matches_traceback
                i = i-1
                j = j-1
            elif traceback_matrix[i][j] == "up":
                seq1_traceback = seq1[i-1] + seq1_traceback
                seq2_traceback = "-" + seq2_traceback
                matches_traceback = " " + matches_traceback
                i = i-1
            elif traceback_matrix[i][j] == "left":
                seq1_traceback = "-" + seq1_traceback
                seq2_traceback = seq2[j-1] + seq2_traceback
                matches_traceback = " " + matches_traceback
                j = j-1

        """ Writes the non-aligned part of the sequences:"""
        # ... before the alignment
        elif i>0:
            if alignment in ["global", "semi"]:
                seq1_traceback = seq1[i-1] + seq1_traceback
                seq2_traceback = "-" + seq2_traceback
                matches_traceback = " " + matches_traceback
                i = i-1
        elif j>0:
            if alignment in ["global", "semi"]:
                seq2_traceback = seq2[j-1] + seq2_traceback
                seq1_traceback = "-" + seq1_traceback
                matches_traceback = " " + matches_traceback
                j = j-1

    # ... after the alignment
    if alignment == "semi":
        spaces_seq1 = len(seq1) - i_max
        spaces_seq2 = len(seq2) - j_max

        for pos in range(1, spaces_seq1):
            seq1_traceback = seq1_traceback + seq1[len(seq1)-pos]
            matches_traceback = matches_traceback + " "
            seq2_traceback = seq2_traceback + "-"
        for pos in range(1, spaces_seq2):
            seq2_traceback = seq2_traceback + seq2[len(seq2)-pos]
            matches_traceback = matches_traceback + " "
            seq1_traceback = seq1_traceback + "-"

    """ Writes the alignment as a string to be printed"""
    score = "score = " + str(score_matrix[i_max+1][j_max+1])

    aligned_sequences = seq1_traceback + "\n" + matches_traceback + \
                        "\n" + seq2_traceback + "\n" + score

    return aligned_sequences

def print_matrix(seq1, seq2, score_matrix):
    """ Prints the score matrix"""

    columns = len(seq2)+2
    rows = len(seq1)+2

    score_matrix_names = [] + score_matrix

    # Add sequence names to the matrix
    score_matrix_names.insert(0, [""] + ["-"] + list(seq2))
    score_matrix_names[1].insert(0, "-")
    for i in range(2,rows):
            score_matrix_names[i].insert(0, seq1[i-2])

    # Converts the list of lists into a viewer friendly matrix
    show_matrix = ""
    for i in range(rows):
        for j in range(columns):
            add_spaces = 5 - len(str(score_matrix_names[i][j]))
            show_matrix = show_matrix + " "*add_spaces + str(score_matrix_names[i][j])
        show_matrix = show_matrix + "\n"

    print show_matrix

    # This fixes a bug, otherwise the score_matrix is changed
    for i in range(0, rows-1):
        score_matrix[i].pop(0)

# main function:
def main():
    # get command line options
    options = parse_commandline()

    # set substitution matrix:
    if options.exchange_matrix == "pam250":
        exchangeMatrix = pam250
    elif options.exchange_matrix == "blosum62":
        exchangeMatrix = blosum62
    elif options.exchange_matrix == "identity":
        exchangeMatrix = identity
    else:
        print "unknown exchange matrix", options.exchange_matrix
        exit(-1)

    # read sequences from fasta file, and catch error reading file
    try:
        sequences = readSequences(open(options.fasta))
    except IOError:
        print "ERROR: cannot open or read fasta input file:", fastafile
        exit(-1)

    # Choose to print score matrix only if -v argument was given
    if options.print_scorematrix:
        printing = True
    else:
        printing = False

    # Prints the sequences in the screen
    for seq in sequences:
        print seq
    print "\n"

    # call alignment routine(s):
    if options.align_global:
        alignment = "global"
        do_global_alignment(sequences, exchangeMatrix, options.gap_penalty, alignment, printing)
    elif options.align_local:
        alignment = "local"
        do_local_alignment(sequences, exchangeMatrix, options.gap_penalty, alignment, printing)
    elif options.align_semiglobal:
        alignment = "semi"
        do_semiglobal_alignment(sequences, exchangeMatrix, options.gap_penalty, alignment, printing)
    else:
        print "BUG! this should not happen."
        exit(-1)


if __name__ == "__main__":
    main()

# last line
