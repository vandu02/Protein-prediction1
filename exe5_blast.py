import numpy as np

from pathlib import Path
from collections import defaultdict, Counter



"""
ATTENTION: Use the following dictionaries to get the correct index for each
           amino acid when accessing any type of matrix (PSSM or substitution
           matrix) parameters. Failure to do so will most likely result in not
           passing the tests.
"""
ALPHABET = 'ACDEFGHIKLMNPQRSTVWY'
AA_TO_INT = {aa: index for index, aa in enumerate(ALPHABET)}
INT_TO_AA = {index: aa for index, aa in enumerate(ALPHABET)}


class BlastDB:

    def __init__(self):
        """
        Initialize the BlastDB class.
        """
        #pass

        self.sequence_indexer = {}
        self.sequence_store = []
        #self.sequence_indexer2 = defaultdict(list)

    def add_sequence(self, sequence):
        """
        Add a sequence to the database.

        :param sequence: a protein sequence (string).
        """
        #pass

        self.sequence_store.append(sequence)
        sequence_index = len(self.sequence_store) - 1
        for i in range(len(sequence) - 2):
            word = sequence[i:i + 3]
            if word in self.sequence_indexer:
                self.sequence_indexer[word].add(sequence_index)
                #self.sequence_indexer2[word].append(sequence_index)
            else:
                self.sequence_indexer[word] = {sequence_index}
                

    def get_sequences(self, word):
        """
        Return all sequences in the database containing a given word.

        :param word: a word (string).

        :return: List with sequences.
        """
        #return ['NOPE']

        if word in self.sequence_indexer:
            return [self.sequence_store[i] for i in self.sequence_indexer[word]]
        else:
            return []


    def get_db_stats(self):
        """
        Return some database statistics:
            - Number of sequences in database
            - Number of different words in database
            - Average number of words per sequence (rounded to nearest int)
            - Average number of sequences per word (rounded to nearest int)

        :return: Tuple with four integer numbers corrsponding to the mentioned
                 statistics (in order of listing above).
        """
        #return (1, 2, 3, 4)

        sequence_count = len(self.sequence_store)
        word_count = len(self.sequence_indexer)
        index_count = sum(len(indices) for indices in self.sequence_indexer.values())

        avg_words_sequence = index_count / sequence_count
        avg_sequences_word = index_count / word_count

        avg_words_sequence = int(round(avg_words_sequence))
        avg_sequences_word = int(round(avg_sequences_word))

        return sequence_count, word_count, avg_words_sequence, avg_sequences_word

class Blast:

    def __init__(self, substitution_matrix):
        """
        Initialize the Blast class with the given substitution_matrix.

        :param substitution_matrix: 20x20 amino acid substitution score matrix.
        """
        #pass

        self.substitution_matrix = substitution_matrix
        self.database = BlastDB()
        self.words = {}
        self.match_dict = {}

        


    def get_words(self, *, sequence=None, pssm=None, T=11):
        """
        Return all words with score >= T for given protein sequence or PSSM.
        Only a sequence or PSSM will be provided, not both at the same time.
        A word may only appear once in the list.

        :param sequence: a protein sequence (string).
        :param pssm: a PSSM (Lx20 matrix, where L is length of sequence).
        :param T: score threshold T for the words.

        :return: List of unique words.
        """
        #return ['AAA', 'YYY']

        words_dict = defaultdict(set)
    
        if sequence and not pssm:
           for pos in range(len(sequence) - 2):
               word = sequence[pos:pos + 3]
               rows = self.substitution_matrix[[AA_TO_INT[aa] for aa in word], :]
               for i, score_i in enumerate(rows[0]):
                   for j, score_j in enumerate(rows[1]):
                       for k, score_k in enumerate(rows[2]):
                           score = score_i + score_j + score_k
                           if score >= T:
                               full_word = f'{INT_TO_AA[i]}{INT_TO_AA[j]}{INT_TO_AA[k]}'
                               #if full_word not in words_dict:
                                #   words_dict[full_word] = set()
                               words_dict[full_word].add((pos, score))

        

        elif not sequence and pssm is not None:
            row_count = pssm.shape[0]
            for pos in range(row_count - 2):
                rows = pssm[pos:pos + 3, :]
                for i, score_i in enumerate(rows[0]):
                    for j, score_j in enumerate(rows[1]):
                        for k, score_k in enumerate(rows[2]):
                            score = score_i + score_j + score_k
                            if score >= T:
                                word = f'{INT_TO_AA[i]}{INT_TO_AA[j]}{INT_TO_AA[k]}'
                                #if word not in words_dict:
                                #    words_dict[word] = set()
                                words_dict[word].add((pos, score))

        self.words = words_dict

        return list(words_dict.keys())

    
    class HSP:
        def __init__(self, sequence='', score=0, query_pos=0, match_pos=0, match='', query='', matrix=[], X=0):
            self.sequence = sequence
            self.score = score
            self.query_pos = query_pos
            self.match_pos = match_pos
            self.match = match
            self.query = query
            self.matrix = matrix
            self.X = X
            
        
        def right_expansion(self):
            self.max_score = self.score
            self.max_sequence=self.sequence
            n = 3
            right_ending = len(self.query) if self.query else self.matrix.shape[0]
            
            while self.query_pos + n < right_ending and self.match_pos + n < len(self.match):
                row_index = AA_TO_INT[self.query[self.query_pos + n]] if self.query else self.query_pos + n #else case if pssm is provided
                current_match = self.match[self.match_pos + n]

                self.score += self.matrix[row_index][AA_TO_INT[current_match]]
                self.sequence = str(self.sequence) + str(current_match)

                if (self.score > self.max_score) or (self.score == self.max_score and len(self.sequence) < len(self.max_sequence)):
                    self.max_score = self.score
                    self.max_sequence = self.sequence

                n += 1
                if self.score <= self.max_score - self.X:
                    break

            self.score = self.max_score
            self.sequence = self.max_sequence

            #return self.sequence, self.score, self.query_pos, self.match_pos, self.match, self.query, self.matrix, self.X
            
        
        def left_expansion(self):
            n = -1
            self.max_score = self.score
            self.max_sequence = self.sequence #initially obtained from right expansion
            max_query_pos = -1
            max_match_pos = -1

            while self.query_pos + n >= 0 and self.match_pos + n >= 0:
                #max_sequence_copy = self.sequence # sequence to keep expanding til we find a better global maximum or drop to S
                max_sequence_ = self.max_sequence # store temporary global maximum
                max_score_ = self.max_score 

                #Find residues to align and their scores in substitution matrix
                row_index = AA_TO_INT[self.query[self.query_pos + n]] if self.query else self.query_pos + n
                current_match = self.match[self.match_pos + n]

                # Extending to the left
                self.score += self.matrix[row_index][AA_TO_INT[current_match]]
                self.sequence = str(current_match) + str(self.sequence)

                #check if score after extension is the new global maximum score, if so, update global values
                if (self.score > max_score_) or (self.score == max_score_ and len(self.sequence) < len(max_sequence_)):
                    max_score_ = self.score
                    max_sequence_ = self.sequence # new global maximum

               # if max_sequence_copy != max_sequence_: #if local and global maximum differ after extension
                    max_query_pos = self.query_pos + n #store initial positions of global maximum in query
                    max_match_pos = self.match_pos + n #store initial positions of global maximum in match
                    self.max_sequence = max_sequence_  # Update new global maximum in self.max_sequence
                    self.max_score = max_score_

                elif self.score <= max_score_ - self.X:
                    break

                n -= 1

            self.sequence = self.max_sequence
            self.score = self.max_score
            if max_query_pos != -1 and max_match_pos != -1:
                self.query_pos, self.match_pos = max_query_pos, max_match_pos



    def find_occurrences(self, word, target):
        i = target.find(word)
        while i != -1:
            yield i 
            i = target.find(word, i+1)


    
    

    def search_one_hit(self, blast_db, *, query=None, pssm=None, T=13, X=5, S=30):
        """
        Search a database for target sequences with a given query sequence or
        PSSM. Return a dictionary where the keys are the target sequences for
        which HSPs have been found and the corresponding values are lists of
        tuples. Each tuple is a HSP with the following elements (and order):
            - Start position of HSP in query sequence
            - Start position of HSP in target sequence
            - Length of the HSP
            - Total score of the HSP
        The same HSP may not appear twice in the list (remove duplictes).
        Only a sequence or PSSM will be provided, not both at the same time.

        :param blast_db: BlastDB class object with protein sequences.
        :param query: query protein sequence.
        :param pssm: query PSSM (Lx20 matrix, where L is length of sequence).
        :param T: score threshold T for the words.
        :param X: drop-off threshold X during extension.
        :param S: score threshold S for the HSP.

        :return: dictionary of target sequences and list of HSP tuples.
        """

        HSPs = defaultdict(set) #Creating a dictionary
        query_words = self.get_words(sequence=query, pssm=pssm, T=T)

        
        if query:
            matrix = self.substitution_matrix
        else:
            matrix = pssm

       
        self.match_dict = self.words.copy() # dictionary of words - matches
        self.matches = self.words.copy() # dictionary of words - positions of every match


        for word in query_words:
            matches = blast_db.get_sequences(word)
            for query_pos, score in self.words[word]:
                for match in matches:
                    match_positions = list(self.find_occurrences(word, match))
                    for match_pos in match_positions:
                        #print('*************** ', match_pos)
                        hsp = self.HSP(sequence=word, query=query, query_pos=query_pos, matrix=matrix, X=X, match=match, match_pos=match_pos, score=score)
                        #print('before right expansion')
                        hsp.right_expansion()
                        #print('right expansion')
                        hsp.left_expansion()
                        #print("VIVA LA MEXICO")
                        #if hsp.match not in HSPs:
                        #        HSPs[hsp.match] = set()
                        if hsp.score >= S:
                            HSPs[hsp.match].add((hsp.query_pos, hsp.match_pos, len(hsp.sequence), hsp.score))

                        #HSPs = {target: hsps for target, hsps in HSPs.items() if hsps}
        

        #print(HSPs)
       
        return HSPs




        """
        import sys
        filename = "output_values.txt"
        with open(filename, "w") as file:
            sys.stdout = file  # Redirect standard output to the file
            print(self.words)        # Print statement will write to the file

        sys.stdout = sys.__stdout__
        """
         
    def get_seq_words(self, sequence):
        sequence_dict = defaultdict(set)
        for n in range(len(sequence) - 2):
            word = sequence[n:n + 3]
            sequence_dict[sequence].add(word)

        return sequence_dict



        

        

    def search_two_hit(self, blast_db, *, query=None, pssm=None, T=11, X=5, S=30, A=40):
        """
        Search a database for target sequences with a given query sequence or
        PSSM. Return a dictionary where the keys are the target sequences for
        which HSPs have been found and the corresponding values are lists of
        tuples. Each tuple is a HSP with the following elements (and order):
            - Start position of HSP in query sequence
            - Start position of HSP in target sequence
            - Length of the HSP
            - Total score of the HSP
        The same HSP may not appear twice in the list (remove duplictes).
        Only a sequence or PSSM will be provided, not both at the same time.

        :param blast_db: BlastDB class object with protein sequences.
        :param query: query protein sequence.
        :param pssm: query PSSM (Lx20 matrix, where L is length of sequence).
        :param T: score threshold T for the words.
        :param X: drop-off threshold X during extension.
        :param S: score threshold S for the HSP.
        :param A: max distance A between two hits for the two-hit method.

        :return: dictionary of target sequences and list of HSP tuples.
        """


        HSPs = defaultdict(set)
        query_words = self.get_words(sequence=query, pssm=pssm, T=T)
        query_set = set(query_words)
        logger = defaultdict(lambda: defaultdict(list))

        if query:
            matrix = self.substitution_matrix
        else:
            matrix = pssm

        match_indices = set().union(*[blast_db.sequence_indexer[word] for word in query_words])

        for n in match_indices:
            match = blast_db.sequence_store[n]
            hit_diagonal = defaultdict(list)
            seq_matches = self.get_seq_words(match)
            hits = set(query_words).intersection(seq_matches[match])

            if len(hits) < 2:
                continue

            pos_seq_words = [(query_pos, score, word) for word in hits for query_pos, score in self.words[word]]
            pos_seq_words.sort(key=lambda a: a[0])

            for query_pos, score, word in pos_seq_words:
                position_hit = self.find_occurrences(word, match)

                for match_pos in position_hit:
                    diagonal = match_pos - query_pos
                    right_hit = (query_pos, match_pos, word, score)
                    hsp_to_add = False
                    remove = -1

                    if hit_diagonal[diagonal]:
                        for a, left_hit in enumerate(hit_diagonal[diagonal]):
                            if hsp_to_add:
                                break

                            if right_hit[0] - left_hit[0] >= 3 and right_hit[1] - left_hit[1] >= 3:
                                distance = right_hit[0] - left_hit[0]

                                if abs(distance) <= A and distance == right_hit[1] - left_hit[1]:
                                    included = any(hsp[1] <= left_hit[1] < right_hit[1] < hsp[1] + hsp[2]
                                                for hsp in logger[match][diagonal])

                                    if not included:
                                        hsp = self.HSP(sequence=word, query=query, query_pos=query_pos, matrix=matrix,
                                                    X=X, match=match, match_pos=match_pos, score=score)
                                        hsp.left_expansion()

                                        if hsp.match_pos <= left_hit[1] + 2:
                                            hsp_query_pos = hsp.query_pos
                                            hsp_match_pos = hsp.match_pos
                                            hsp.query_pos = right_hit[0]
                                            hsp.match_pos = right_hit[1]
                                            hsp.right_expansion()
                                            logger[match][diagonal].append((hsp_query_pos, hsp_match_pos, len(hsp.sequence)))

                                            if hsp.score >= S:
                                                HSPs[match].add((hsp_query_pos, hsp_match_pos, len(hsp.sequence), hsp.score))
                                                remove = a
                                                hsp_to_add = True

                    if not hsp_to_add:
                        hit_diagonal[diagonal].append(right_hit)
                    else:
                        del hit_diagonal[diagonal][remove]

        return HSPs
