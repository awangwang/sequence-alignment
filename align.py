
"""
Author: Andy Wang

This python file implements a dymamic programming algorithm and recursive traceback
for pseudo-global alignment and local alignment (Smith-Waterman).
Usage: python align.py input_file output_file

"""


import sys
import re


#### ------ USEFUL FUNCTIONS ------- ####
def fuzzy_equals(a, b):
    """
    Checks if two floating point numbers are equivalent.
    """
    epsilon = 10**(-6) 
    return (abs(a - b) < epsilon)


#### ------- CLASSES ------- ####

class MatchMatrix(object):
    """
    Match matrix class stores the scores of matches in a data structure
    """

    def __init__(self):
        self.match_matrix = {}

    def set_score(self, a, b, score):
        """
        Updates or adds a score for a specified match
        Uses tuple of a and b as key to store the value of score
        Tuples are ordered, so this data structure can handle asymmetric match matrices

        Inputs:
           a = the character from sequence A
           b = the character from sequence B
           score = the score to set it for
        """

        self.match_matrix[(a, b)] = score

    def get_score(self, a, b):
        """
        Returns the score for a particular match, where a is the
        character from sequence a and b is from sequence b.

        Inputs:
           a = the character from sequence A
           b = the character from sequence B
        Returns:
           the score of that match
        """

        return self.match_matrix[(a, b)]


class ScoreMatrix(object):
    """
    Object to store a score matrix, which generated during the alignment process. The score matrix consists of a 2-D array of
    ScoreEntries that are updated during alignment and used to output the maximum alignment.
    """

    def __init__(self, name, nrow, ncol):
        self.name = name # identifier for the score matrix - Ix, Iy, or M
        self.nrow = nrow
        self.ncol = ncol
        self.score_matrix = [[{} for _ in range(self.ncol)] for _ in range(self.nrow)]

    def initialize_matrix(self, nrow, ncol):
        """
        (Re)initializes the matrix once all the alignment parameters are loaded
        in from the input file.

        Inputs:
           nrow = number of rows in the matrix
           ncol = number of columns in the matrix
        """

        self.nrow = nrow
        self.ncol = ncol
        self.score_matrix = [[{} for _ in range(self.ncol)] for _ in range(self.nrow)]

        # set scores of 0th row and column to 0
        for i in range(self.nrow):
            self.set_score(i, 0, 0)
        for j in range(self.ncol):
            self.set_score(0, j, 0)

    def get_score(self, row, col):
        """
        Returns the score of a cell in the matrix

        Inputs:
           row = index of the row of the cell
           col = index of the column of the cell
        Returns:
            the score of that cell
        """

        return self.score_matrix[row][col]['score']
        
    def set_score(self, row, col, score): 
        """
        Updates or adds a score for a cell in the matrix

        Inputs:
           row = index of the row of the cell
           col = index of the column of the cell
           score = the score to set it for
        """

        self.score_matrix[row][col]['score'] = score

    def get_pointers(self, row, col):
        """
        Returns the indices of the entries that are pointed to
        This should be formatted as a list of tuples:
         ex. [(1,1), (1,0)]
        """
        grid_object = self.score_matrix[row][col]

        # edge case for empty object
        if 'matrix' not in grid_object:
            return []

        pointers = []
        for i in range(len(grid_object['matrix'])):
            pointer = (grid_object['matrix'][i], grid_object['row'][i], grid_object['col'][i])
            pointers.append(pointer)
        return pointers

    def set_pointers(self, row, col, matrix, pointer_row, pointer_col):
        """
        Updates or adds a pointer for a cell in the matrix

        Inputs:
           row = index of the row of the cell
           col = index of the column of the cell
           matrix = matrix the pointer is pointing to
           pointer_row = index of the row the pointer is pointing to
           pointer_col = index of the column the pointer is pointing to
        """

        grid_object = self.score_matrix[row][col]

        if 'matrix' not in grid_object:
            grid_object['matrix'] = [matrix]
            grid_object['row'] = [pointer_row]
            grid_object['col'] = [pointer_col]
        else:
            grid_object['matrix'].append(matrix)
            grid_object['row'].append(pointer_row)
            grid_object['col'].append(pointer_col)

    def print_scores(self):
        """
        Returns a nicely formatted string containing the scores in the score matrix. Use this for debugging!

        Example:
        M=
            0.0, 0.0, 0.0, 0.0, 0.0
            0.0, 1.0, 0.0, 0.0, 0.0
            0.0, 1.0, 1.0, 1.0, 1.0
            0.0, 0.0, 1.0, 1.0, 1.0
            0.0, 0.0, 2.0, 2.0, 1.0
            0.0, 0.0, 1.0, 2.0, 3.0

        """

        scores = ""
        for i in range(self.nrow):
            for j in range(self.ncol):
                if j != 0:
                    scores += ", "
                scores += str(round(self.score_matrix[i][j]['score'], 1))
            scores += "\n"
        return scores

    def print_pointers(self):
        """
        Returns a nicely formatted string containing the pointers for each entry in the score matrix. Use this for debugging!
        """

        pointers = ""
        for i in range(1, self.nrow):
            for j in range(1, self.ncol):
                matrix = self.score_matrix[i][j]['matrix']
                row = self.score_matrix[i][j]['row']
                col = self.score_matrix[i][j]['col']
                pointers += f"({matrix}, {row}, {col})"
            pointers += "\n"
        return pointers


class AlignmentParameters(object):
    """
    Object to hold a set of alignment parameters from an input file.
    """

    def __init__(self):
        # default values for variables that are filled in by reading
        # the input alignment file
        self.seq_a = ""
        self.seq_b = ""
        self.len_seq_a = 0
        self.len_seq_b = 0
        self.global_alignment = False 
        self.dx = 0
        self.ex = 0
        self.dy = 0
        self.ey = 0
        self.alphabet_a = "" 
        self.alphabet_b = ""
        self.len_alphabet_a = 0
        self.len_alphabet_b = 0
        self.match_matrix = MatchMatrix()

    def load_params_from_file(self, input_file): 
        """
        Reads the parameters from an input file and stores in the object

        Input:
           input_file = specially formatted alignment input file
        """

        with open(input_file, 'r') as file:
            lines = file.readlines()
        
        self.seq_a = lines[0].strip()
        self.seq_b = lines[1].strip()
        self.len_seq_a = len(self.seq_a)
        self.len_seq_b = len(self.seq_b)
        self.global_alignment = (int(lines[2]) == 0)
        penalties = [float(penalty) for penalty in lines[3].split()]
        self.dx = penalties[0]
        self.ex = penalties[1]
        self.dy = penalties[2]
        self.ey = penalties[3]
        self.len_alphabet_a = int(lines[4])
        self.alphabet_a = lines[5].strip()
        self.len_alphabet_b = int(lines[6])
        self.alphabet_b = lines[7].strip()

        for line in lines[8:]:
            cols = line.split()
            if len(cols) < 5:
                break
            a = cols[2]
            b = cols[3]
            score = float(cols[4])
            self.match_matrix.set_score(a, b, score)


class Align(object):
    """
    Object to hold and run an alignment; running is accomplished by using "align()"
    """

    def __init__(self, input_file, output_file):
        """
        Input:
            input_file = file with the input for running an alignment
            output_file = file to write the output alignments to
        """
        self.input_file = input_file
        self.output_file = output_file
        self.align_params = AlignmentParameters() 

        self.m_matrix = ScoreMatrix('M', 0, 0)
        self.ix_matrix = ScoreMatrix('Ix', 0, 0)
        self.iy_matrix = ScoreMatrix('Iy', 0, 0)
        self.matrix_dictionary = {'M': self.m_matrix, 'Ix': self.ix_matrix, 'Iy': self.iy_matrix}

    def align(self):
        """
        Main method for running alignment.
        """

        # load the alignment parameters into the align_params object
        self.align_params.load_params_from_file(self.input_file)

        # populate the score matrices based on the input parameters
        self.populate_score_matrices()

        # perform a traceback and write the output to an output file
        max_val, max_loc = self.find_traceback_start()
        all_alignments = self.traceback(max_loc)
        self.write_output(max_val, all_alignments)

    def populate_score_matrices(self):
        """
        Method to populate the score matrices based on the data in align_params.
        Should call update(i,j) for each entry in the score matrices
        """

        self.m_matrix.initialize_matrix(self.align_params.len_seq_a + 1, self.align_params.len_seq_b + 1)
        self.ix_matrix.initialize_matrix(self.align_params.len_seq_a + 1, self.align_params.len_seq_b + 1)
        self.iy_matrix.initialize_matrix(self.align_params.len_seq_a + 1, self.align_params.len_seq_b + 1)

        for i in range(1, self.align_params.len_seq_a + 1):
            for j in range(1, self.align_params.len_seq_b + 1):
                self.update(i, j)

    def update(self, row, col):
        """
        Method to update the matrices at a given row and column index.

        Input:
           row = the row index to update
           col = the column index to update
        """

        self.update_m(row, col)
        self.update_ix(row, col)
        self.update_iy(row, col)

    def update_m(self, row, col):
        """
        Updates score and pointers of M matrix at a given row and column index

        Input:
           row = the row index to update
           col = the column index to update
        """

        score_i_j = self.align_params.match_matrix.get_score(self.align_params.seq_a[row - 1], self.align_params.seq_b[col - 1])
        m_candidate = self.m_matrix.get_score(row - 1, col - 1) + score_i_j
        ix_candidate = self.ix_matrix.get_score(row - 1, col - 1) + score_i_j
        iy_candidate = self.iy_matrix.get_score(row - 1, col - 1) + score_i_j
        max_score = max(m_candidate, ix_candidate, iy_candidate)

        if self.align_params.global_alignment == False and max_score < 0:
            self.m_matrix.set_score(row, col, 0)
        else:
            self.m_matrix.set_score(row, col, max_score)

        if fuzzy_equals(m_candidate, max_score):
            self.m_matrix.set_pointers(row, col, 'M', row - 1, col - 1)
        if fuzzy_equals(ix_candidate, max_score):
            self.m_matrix.set_pointers(row, col, 'Ix', row - 1, col - 1)
        if fuzzy_equals(iy_candidate, max_score):
            self.m_matrix.set_pointers(row, col, 'Iy', row - 1, col - 1)

    def update_ix(self, row, col):
        """
        Updates score and pointers of Ix matrix at a given row and column index

        Input:
           row = the row index to update
           col = the column index to update
        """

        m_candidate = self.m_matrix.get_score(row - 1, col) - self.align_params.dy
        ix_candidate = self.ix_matrix.get_score(row - 1, col) - self.align_params.ey
        max_score = max(m_candidate, ix_candidate)

        if self.align_params.global_alignment == False and max_score < 0:
            self.ix_matrix.set_score(row, col, 0)
        else:
            self.ix_matrix.set_score(row, col, max_score)

        if fuzzy_equals(m_candidate, max_score):
            self.ix_matrix.set_pointers(row, col, 'M', row - 1, col)
        if fuzzy_equals(ix_candidate, max_score):
            self.ix_matrix.set_pointers(row, col, 'Ix', row - 1, col)

    def update_iy(self, row, col):
        """
        Updates score and pointers of Iy matrix at a given row and column index

        Input:
           row = the row index to update
           col = the column index to update
        """

        m_candidate = self.m_matrix.get_score(row, col - 1) - self.align_params.dx
        iy_candidate = self.iy_matrix.get_score(row, col - 1) - self.align_params.ex
        max_score = max(m_candidate, iy_candidate)

        if self.align_params.global_alignment == False and max_score < 0:
            self.iy_matrix.set_score(row, col, 0)
        else:
            self.iy_matrix.set_score(row, col, max_score)

        if fuzzy_equals(m_candidate, max_score):
            self.iy_matrix.set_pointers(row, col, 'M', row, col - 1)
        if fuzzy_equals(iy_candidate, max_score):
            self.iy_matrix.set_pointers(row, col, 'Iy', row, col - 1)

    def find_traceback_start(self):
        """
        Finds the location to start the traceback..
        Think carefully about how to set this up for local 

        Returns:
            (max_val, max_loc) where max_val is the best score
            max_loc is a set() containing tuples with the (i,j) location(s) to start the traceback
             (ex. [(1,2), (3,4)])
        """

        nrow = self.align_params.len_seq_a + 1
        ncol = self.align_params.len_seq_b + 1

        max_val = -sys.float_info.max
        max_loc = set()

        # pseudo-global: find highest score in last row or column of M
        if self.align_params.global_alignment:
            for i in range(nrow):
                score = self.m_matrix.get_score(i, ncol - 1)
                if score > max_val:
                    max_val = score
                    max_loc = set()
                    max_loc.add(('M', i, ncol - 1))
                elif fuzzy_equals(score, max_val):
                    max_loc.add(('M', i, ncol - 1))
            
            for j in range(ncol):
                score = self.m_matrix.get_score(nrow - 1, j)
                if score > max_val:
                    max_val = score
                    max_loc = set()
                    max_loc.add(('M', nrow - 1, j))
                elif fuzzy_equals(score, max_val):
                    max_loc.add(('M', nrow - 1, j))
        
        # local: find highest score anywhere in M
        else:
            for i in range(1, nrow):
                for j in range(1, ncol):
                    score = self.m_matrix.get_score(i, j)
                    if score > max_val:
                        max_val = score
                        max_loc = set()
                        max_loc.add(('M', i, j))
                    elif fuzzy_equals(score, max_val):
                        max_loc.add(('M', i, j))

        return max_val, max_loc

    def traceback(self, max_loc):
        """
        Performs a traceback

        Input:
            max_loc = starting location(s) with highest score
        """

        all_alignments = []
        for start_loc in max_loc:
            alignments_at_loc = set(self.rec_trace(start_loc[0], start_loc[1], start_loc[2]))
            all_alignments.extend(alignments_at_loc)
        trimmed_alignments = self.trim_end_gaps(all_alignments)
        return trimmed_alignments

    def rec_trace(self, matrix_name, row, col):
        """
        Recursive helper for traceback()

        Input:
            matrix_name = name of the matrix that the traceback is currently in
            row = row that the traceback is currently in
            col = column that the traceback is currently in
        """

        matrix = self.matrix_dictionary[matrix_name]
        
        # base cases:
        # terminate if the 0th row or column is reached
        if row == 0:
            return [("", "")]
        if col == 0:
            return [("", "")]
        # local alignment: terminate at cells with score of 0
        if not self.align_params.global_alignment and matrix.get_score(row, col) == 0:
            return [("", "")]
        
        emits_to_return = []
        pointers = matrix.get_pointers(row, col)

        for pointer in pointers:
            # emit both letters
            if (pointer[1] == row - 1) and (pointer[2] == col - 1):
                emit_a = self.align_params.seq_a[row - 1]
                emit_b = self.align_params.seq_b[col - 1]

            # gap in A
            if (pointer[1] == row) and (pointer[2] == col - 1):
                emit_a = '_'
                emit_b = self.align_params.seq_b[col - 1]

            # gap in B
            if (pointer[1] == row - 1) and (pointer[2] == col):
                emit_a = self.align_params.seq_a[row - 1]
                emit_b = '_'
            
            recursed_emits = self.rec_trace(pointer[0], pointer[1], pointer[2])
            for subseq in recursed_emits:
                emits_to_return.append((emit_a + subseq[0], emit_b + subseq[1]))
        
        return emits_to_return

    def trim_end_gaps(self, all_alignments):
        """
        Removes end gaps from aligned sequences

        Input:
            all_alignments = list of all optimal sequence alignments
        Returns:
            the aligned sequences with end gaps removed
        """

        trimmed_alignments = []

        for alignment in all_alignments:
            start_gaps = max(self.count_gaps(alignment[0], 'start'), self.count_gaps(alignment[1], 'start'))
            end_gaps = max(self.count_gaps(alignment[0], 'end'), self.count_gaps(alignment[1], 'end'))
            
            if end_gaps == 0:
                trimmed_alignment = (alignment[0][start_gaps:], alignment[1][start_gaps:])
                trimmed_alignments.append(trimmed_alignment)
            else:
                trimmed_alignment = (alignment[0][start_gaps: -end_gaps], alignment[1][start_gaps: -end_gaps])
                trimmed_alignments.append(trimmed_alignment)
        
        return trimmed_alignments

    def count_gaps(self, seq, direction):
        """
        Counts the number of gaps/underscores at either end of a sequence, depending on
        the value of `direction`, using regex

        Input:
            seq = the sequence being counted
            direction = either 'start' or 'end', indicating the direction to count from
        Returns:
            the number of gaps at the specified side of seq
        """

        if direction == 'start':
            gaps = re.search(r'^_+', seq)
            if gaps:
                return len(gaps.group(0))
            return 0
        elif direction == 'end':
            gaps = re.search(r'_+$', seq)
            if gaps:
                return len(gaps.group(0))
            return 0

    def write_output(self, max_val, all_alignments):
        """
        Writes the best score as well as all the sequence alignments that have that score
        to the specified output file

        Input:
            max_val = the best/winning alignment score
            all_alignments = list of all optimal sequence alignments
        """

        with open(self.output_file, 'w') as file:
            file.write(str(round(max_val, 1)) + '\n')
            for alignment in all_alignments:
                file.write('\n')
                alignment_a = alignment[0][::-1]
                alignment_b = alignment[1][::-1]
                file.write(alignment_a + '\n')
                file.write(alignment_b + '\n')


def main():
    # check that the file is being properly used
    if (len(sys.argv) !=3):
        print("Please specify an input file and an output file as args.")
        return

    # input variables
    input_file = sys.argv[1]
    output_file = sys.argv[2]

    # create an align object and run
    align = Align(input_file, output_file)
    align.align()


if __name__=="__main__":
    main()
