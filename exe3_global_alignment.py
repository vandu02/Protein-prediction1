import numpy as np

class GlobalAlignment:
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
        self.substituion_matrix = matrix
        self.score_matrix = np.zeros((len(string2) + 1, len(string1) + 1), dtype=np.int32)
        self.align()

    def align(self):
        """
        Align given strings using the Needleman-Wunsch algorithm,
        store the alignments and the score matrix used to compute those alignments.
        NB: score matrix and the substitution matrix are different matrices!
        """

        for i in range(len(self.string2) + 1):
            self.score_matrix[i][0] = self.gap_penalty * i

        for j in range(len(self.string1) + 1):
            self.score_matrix[0][j] = self.gap_penalty * j

        for i in range(1, len(self.string2) + 1):
            for j in range(1, len(self.string1) + 1):
                diagonal_score = self.score_matrix[i - 1][j - 1] + self.substituion_matrix[self.string2[i - 1]][self.string1[j - 1]]
                up_score = self.score_matrix[i - 1][j] + self.gap_penalty
                left_score = self.score_matrix[i][j - 1] + self.gap_penalty

                self.score_matrix[i][j] = max(diagonal_score, up_score, left_score)


    def get_best_score(self):
        """
        :return: the highest score for the aligned strings, int

        """
        #return 4
        best_score = self.score_matrix[-1][-1]
        return best_score

    def get_number_of_alignments(self):
        """
        :return: number of found alignments with the best score
        """
        #return 43

        best_score = self.get_best_score()
        count = 0

        def backtrack(i, j):
            nonlocal count

            if i == 0 and j == 0:
                count += 1
                return

            if i > 0 and j > 0 and self.score_matrix[i][j] == self.score_matrix[i - 1][j - 1] + self.substituion_matrix[self.string2[i - 1]][self.string1[j - 1]]:
                backtrack(i - 1, j - 1)

            if i > 0 and self.score_matrix[i][j] == self.score_matrix[i - 1][j] + self.gap_penalty:
                backtrack(i - 1, j)

            if j > 0 and self.score_matrix[i][j] == self.score_matrix[i][j - 1] + self.gap_penalty:
                backtrack(i, j - 1)

        backtrack(len(self.string2), len(self.string1))
        return count


    def get_alignments(self):
        """
        :return: list of alignments, where each alignment is represented
                 as a tuple of aligned strings
        """
        #return [
        #    ('ADMI-NS', 'ADMIRES'), ('ADMIN-S', 'ADMIRES')
        #]
        
        best_score = self.get_best_score()
        alignments = []

        def backtrack(i, j, align1, align2):
            if i == 0 and j == 0:
                alignments.append((align1[::-1], align2[::-1]))
                return

            if i > 0 and j > 0 and self.score_matrix[i][j] == self.score_matrix[i - 1][j - 1] + self.substituion_matrix[self.string2[i - 1]][self.string1[j - 1]]:
                backtrack(i - 1, j - 1, align1 + self.string1[j - 1], align2 + self.string2[i - 1])

            if i > 0 and self.score_matrix[i][j] == self.score_matrix[i - 1][j] + self.gap_penalty:
                backtrack(i - 1, j, align1 + '-', align2 + self.string2[i - 1])

            if j > 0 and self.score_matrix[i][j] == self.score_matrix[i][j - 1] + self.gap_penalty:
                backtrack(i, j - 1, align1 + self.string1[j - 1], align2 + '-')

        backtrack(len(self.string2), len(self.string1), '', '')

        return alignments

    def get_score_matrix(self):
        """
        :return: matrix built during the alignment process as a list of lists
        """
        """
        return [
            [0, -1, -2, -3, -4, -5, -6],
            [-1, 1, 0, -1, -2, -3, -4],
            [-2, 0, 2, 1, 0, -1, -2],
            [-3, -1, 1, 3, 2, 1, 0],
            [-4, -2, 0, 2, 4, 3, 2],
            [-5, -3, -1, 1, 3, 4, 3],
            [-6, -4, -2, 0, 2, 3, 4],
            [-7, -5, -3, -1, 1, 2, 4]
        ]
        """
        return self.score_matrix.tolist()
