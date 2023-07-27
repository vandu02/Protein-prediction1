import numpy as np

"""
ATTENTION: Use the following dictionaries to get the correct index for each
           amino acid when accessing any type of matrix or array provided as
           parameters. Further, use those indices when generating or returning
           any matrices or arrays. Failure to do so will most likely result in
           not passing the tests.
EXAMPLE: To access the substitution frequency from alanine 'A' to proline 'P'
         in the bg_matrix use bg_matrix[AA_TO_INT['A'], AA_TO_INT['P']].
"""
ALPHABET = 'ACDEFGHIKLMNPQRSTVWY-'
AA_TO_INT = {aa: index for index, aa in enumerate(ALPHABET)}
INT_TO_AA = {index: aa for index, aa in enumerate(ALPHABET)}
GAP_INDEX = AA_TO_INT['-']


class MSA:

    def __init__(self, sequences):
        """
        Initialize the MSA class with the provided list of sequences. Check the
        sequences for correctness. Pre-calculate any statistics you seem fit.

        :param sequences: List containing the MSA sequences.
        """
        #pass

        #Check if sequences are right
        if len(sequences) == 0:
            raise TypeError("No sequences provided.")

        for i, seq in enumerate(sequences):
            if len(seq) == 0:
                raise TypeError(f"Empty sequence found at index {i}.")

            for a in seq:
                if a not in ALPHABET:
                    raise TypeError(f"Invalid character '{a}' in sequence at index {i}.")


        lengths = [len(seq) for seq in sequences]
        equal_lengths = all(length == lengths[0] for length in lengths)
 
        if equal_lengths == False:
            raise TypeError("Sequences have different lenghts.")

        self.sequence = sequences
        self.MSA_len = len(self.sequence[0])
        self.MSA_array = np.array([list(seq) for seq in sequences])



    def get_pssm(self, *, bg_matrix=None, beta=10, use_sequence_weights=False,
                 redistribute_gaps=False, add_pseudocounts=False):
        """
        Return a PSSM for the underlying MSA. Use the appropriate refinements 
        according to the parameters. If no bg_matrix is specified, use uniform 
        background and pair frequencies.
        Every row in the resulting PSSM corresponds to a non-gap position in 
        the primary sequence of the MSA (i.e. the first one).
        Every column in the PSSM corresponds to one of the 20 amino acids.
        Values that would be -inf must be replaced by -20 in the final PSSM.
        Before casting to dtype=numpy.int64, round all values to the nearest
        integer (do not just FLOOR all values).

        :param bg_matrix: Amino acid pair frequencies as numpy array (20, 20).
                          Access the matrix using the indices from AA_TO_INT.
        :param beta: Beta value (float) used to weight the pseudocounts 
                     against the observed amino acids in the MSA.
        :param use_sequence_weights: Calculate and apply sequence weights.
        :param redistribute_gaps: Redistribute the gaps according to the 
                                  background frequencies.
        :param add_pseudocounts: Calculate and add pseudocounts according 
                                 to the background frequencies.

        :return: PSSM as numpy array of shape (L x 20, dtype=numpy.int64).
                 L = ungapped length of the primary sequence.
        """
        #pssm = np.zeros((20, 20))

        #return np.rint(pssm).astype(np.int64)

        pssm = np.zeros((self.MSA_len, 21), dtype=np.float64)
        primary_seq = self.MSA_array[0]

        if bg_matrix is None:
            bg_matrix = np.full((20, 20), 1.0 / (20 * 20)) #uniform distribution
        
        if not use_sequence_weights:
            for i in range(self.MSA_len):
                unique_AAs, counts = np.unique(self.MSA_array[:, i], return_counts=True)
                for j, aa in enumerate(unique_AAs):
                    if aa in AA_TO_INT:
                        pssm[i, AA_TO_INT[aa]] = counts[j]
        else:
            '''
            1) Obtain sequence weights (W0, W1 corresponds to weigths for sequence 0, 1)
            2) While traversing along all alignments at position i, we gather the unique amino acids at this position
                and on which sequence where they found
            '''
            seq_weights = self.get_sequence_weights()
            for i in range(self.MSA_len):
                unique_AAs, seq_indices = np.unique(self.MSA_array[:, i], return_inverse=True)
                for j, aa in enumerate(unique_AAs):
                    if aa in AA_TO_INT:
                        weight_indices = np.where(seq_indices == j)[0]
                        pssm[i, AA_TO_INT[aa]] = np.sum(seq_weights[weight_indices])

        if redistribute_gaps:
            expected_bg_freq = np.sum(bg_matrix, axis=0)
            pssm[:, :20] = pssm[:, :20] + pssm[:, 20].reshape(-1, 1) * expected_bg_freq.reshape(1, -1)
        
        pssm = pssm[:, :20]
        
        if add_pseudocounts:
            alpha = self.get_number_of_observations() - 1
            pseudocount_mat = np.zeros(pssm.shape)
            p = np.sum(bg_matrix, axis=0)
            for i in range(pssm.shape[0]):
                for j in range(pssm.shape[1]):
                    if pssm[i, j] != 0.0:
                        for a in range(20):
                            pseudocount_mat[i, a] += bg_matrix[j, a] * (pssm[i, j] / p[j])
            pssm = (alpha * pssm + float(beta) * pseudocount_mat) / (alpha + beta)

        row_sums = pssm.sum(axis=1)
        pssm = pssm / row_sums[:, np.newaxis]

        p = np.sum(bg_matrix, axis=0)
        pssm = pssm / p.reshape((1, -1))

        pssm = 2 * np.log2(pssm + 1e-20)

        pssm[pssm < -20] = -20
        
        target_row_idxes = np.where(primary_seq == '-')[0]
        pssm = np.delete(pssm, target_row_idxes, axis=0)

        pssm = np.rint(pssm).astype(np.int64)
        return pssm


    def get_size(self):
        """
        Return the number of sequences in the MSA and the MSA length, i.e.
        the number of columns in the MSA. This includes gaps.

        :return: Tuple of two integers. First element is the number of
                 sequences in the MSA, second element is the MSA length.
        """
        #return (-1, -1)
        return len(self.sequence), len(self.sequence[0])

    def get_primary_sequence(self):
        """
        Return the primary sequence of the MSA. In this exercise, the primary
        sequence is always the first sequence of the MSA. The returned 
        sequence must NOT include gap characters.

        :return: String containing the ungapped primary sequence.
        """
        #return 'NOPE'

        primary_seq = self.sequence[0].replace('-','')
        return primary_seq

    def get_sequence_weights(self):
        """
        Return the calculated sequence weights for all sequences in the MSA.
        The order of weights in the array must be equal to the order of the
        sequences in the MSA.

        :return: Numpy array (dtype=numpy.float64) containing the weights for
                 all sequences in the MSA.
        """
        #weights = np.zeros(10)

        #return weights.astype(np.float64)

        num_sequences, msa_len = self.MSA_array.shape
        weight_matrix = np.zeros((num_sequences, msa_len), dtype=np.float64)
        unique_counts = np.zeros(msa_len)

        for i in range(msa_len):
            unique_elements, indices = np.unique(self.MSA_array[:, i], return_inverse=True)
            unique_counts[i] = len(unique_elements)

            if unique_counts[i] > 1:
                for u in range(len(unique_elements)):
                    element_indices = np.where(indices == u)
                    element_count = len(element_indices[0])
                    weight_matrix[element_indices, i] = 1.0 / (unique_counts[i] * element_count)

        weights = np.sum(weight_matrix, axis=1) #sum horizontally
        return weights.astype(np.float64)


    def get_number_of_observations(self):
        """
        Return the estimated number of independent observations in the MSA.

        :return: Estimate of independent observation (dtype=numpy.float64).
        """
        #num_obs = np.int64(-1)

        #return num_obs.astype(np.float64)

        num_obs = np.mean([len(np.unique(self.MSA_array[:, i])) for i in range(self.MSA_len)])
        return num_obs.astype(np.float64)
