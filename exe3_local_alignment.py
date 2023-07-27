import numpy as np


class LocalAlignment:
    def __init__(self, string1, string2, gap_penalty, matrix):
        """
        :param string1: first string to be aligned, string
        :param string2: second string to be aligned, string
        :param gap_penalty: gap penalty, integer
        :param matrix: substitution matrix containing scores for amino acid
                       matches and mismatches, dict

        Attention! string1 is used to index columns, string2 is used to index rows
        """
        self.string1 = string1
        self.string2 = string2
        self.gap_penalty = gap_penalty
        self.substitution_matrix = matrix
        self.score_matrix = np.zeros((len(string2) + 1, len(string1) + 1), dtype=np.int32)
        self.align()

    def align(self):
        """
        Align given strings using the Smith-Waterman algorithm.
        NB: score matrix and the substitution matrix are different matrices!
        """

        for i in range(1, len(self.string2) + 1):
            for j in range(1, len(self.string1) + 1):
                match_score = self.substitution_matrix[self.string2[i-1]][self.string1[j-1]]
                diagonal_score = self.score_matrix[i-1][j-1] + match_score
                top_score = self.score_matrix[i-1][j] + self.gap_penalty
                left_score = self.score_matrix[i][j-1] + self.gap_penalty
                self.score_matrix[i][j] = max(0, diagonal_score, top_score, left_score)


    def has_alignment(self):
        """
        :return: True if a local alignment has been found, False otherwise
        """
        #return True
        return np.any(self.score_matrix > 0)

    def get_alignment(self):
        """
        :return: alignment represented as a tuple of aligned strings
        """
        #return ('DAC', 'DAC')
        max_score = np.max(self.score_matrix)
        max_indices = np.where(self.score_matrix == max_score)
        max_row, max_col = max_indices[0][0], max_indices[1][0]

        aligned_string1 = ""
        aligned_string2 = ""
        while max_row > 0 and max_col > 0 and self.score_matrix[max_row][max_col] > 0:
            current_score = self.score_matrix[max_row][max_col]
            diagonal_score = self.score_matrix[max_row - 1][max_col - 1]
            top_score = self.score_matrix[max_row - 1][max_col]
            left_score = self.score_matrix[max_row][max_col - 1]

            if current_score == diagonal_score + self.substitution_matrix[self.string2[max_row-1]][self.string1[max_col-1]]:
                aligned_string1 = self.string1[max_col-1] + aligned_string1
                aligned_string2 = self.string2[max_row-1] + aligned_string2
                max_row -= 1
                max_col -= 1
            elif current_score == top_score + self.gap_penalty:
                aligned_string1 = '-' + aligned_string1
                aligned_string2 = self.string2[max_row-1] + aligned_string2
                max_row -= 1
            elif current_score == left_score + self.gap_penalty:
                aligned_string1 = self.string1[max_col-1] + aligned_string1
                aligned_string2 = '-' + aligned_string2
                max_col -= 1
            else:
                break 

        if not aligned_string1 or not aligned_string2:
            return "", ""

        return aligned_string1, aligned_string2

    def is_residue_aligned(self, string_number, residue_index):
        """
        :param string_number: number of the string (1 for string1, 2 for string2) to check
        :param residue_index: index of the residue to check
        :return: True if the residue with a given index in a given string has been aligned
                 False otherwise
        """
        #return False

        aligned_string1, aligned_string2 = self.get_alignment()
        
        if not aligned_string1 or not aligned_string2:
            return False
        
        #remove underscores from alignments
        align1 = aligned_string1.replace('-', '')
        align2 = aligned_string2.replace('-', '')

        #2 cases
        # str.partition(sep) - string separator!
        # Partition returns a tuple with 3 elements (a,b,c)
        # a = all elements in string1 that are before align1
        # b = align1
        # c = all elements in string1 that are after align1

        '''
        if string_number == 1:
            seq_seps = self.string1.partition(align1)
        else:
            seq_seps = self.string2.partition(align2)
        
        return len(seq_seps[0]) - 1 < residue_index < len(seq_seps[0]) + len(seq_seps[1])
        '''

        
        if string_number == 1:
            seq_seps = self.string1.partition(align1)
            if residue_index >= 0 and residue_index < len(seq_seps[0]) + len(seq_seps[1]):
                return self.string1[residue_index] in aligned_string1
            else:
                return False
            


        else: #string number == 2
            seq_seps = self.string2.partition(align2)       
            if residue_index >= 0 and residue_index < len(seq_seps[0]) + len(seq_seps[1]):
                return self.string2[residue_index] in aligned_string2
            else:
                return False
       
    
        



        """

        if string_number == 1:
            if residue_index >= 0 and residue_index < len(self.string1):
                return self.string1[residue_index] in aligned_string1
            else:
                return False
        elif string_number == 2:
            if residue_index >= 0 and residue_index < len(self.string2):
                return self.string2[residue_index] in aligned_string2
            else:
                return False
        else:
            return False
        """



