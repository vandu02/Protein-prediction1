class Genome:

    def __init__(self, genome):
        """
        Initialize the Genome class with the provided genome sequence.

        :param genome: String with the genome sequence.
        """
        self.genome = genome
        self.length = len(genome)
        self.a_count = genome.count('A')
        self.t_count = genome.count('T') 
        self.g_count = genome.count('G')
        self.c_count = genome.count('C')  
        self.nucleotides = 'ATGC'
        pass

    def get_at_content(self):
        """
        Return the AT content of the genome sequence, i.e. the combined
        fraction of 'A' and 'T' in the entire genome sequence.

        :return: AT content (float, rounded to 6 digits)
        """
        #a_count = self.genome.count('A')
        #t_count = self.genome.count('T')
        
        at_cont = (self.a_count + self.t_count) / self.length
        return round(at_cont, 6)

    def get_codon_dist(self):
        """
        Return the expected codon distribution (fractions) based on the
        distribution (fractions) of the four different nucleotides (ATGC).

        :return: Tree-like structure made out of nested dictionaries. The nodes
                 represent the different nucleotides and the path down the tree
                 forms the corresponding codons. The leafs contain the expected
                 codon frequencies (rounded to 6 digits).
        """
        self.codon_dist = {}
    
        prob_a = self.a_count / self.length
        prob_t = self.t_count / self.length
        prob_g = self.g_count / self.length
        prob_c = self.c_count / self.length
    
        for nucleotide_1 in self.nucleotides:
            self.codon_dist[nucleotide_1] = {}
        
            if nucleotide_1 == "A":
                nuc_freq_1 = prob_a
            elif nucleotide_1 == "T":
                nuc_freq_1 = prob_t
            elif nucleotide_1 == "G":
                nuc_freq_1 = prob_g
            elif nucleotide_1 == "C":
                nuc_freq_1 = prob_c
            
            for nucleotide_2 in self.nucleotides:
                self.codon_dist[nucleotide_1][nucleotide_2] = {}
                
                if nucleotide_2 == "A":
                    nuc_freq_2 = prob_a
                elif nucleotide_2 == "T":
                    nuc_freq_2 = prob_t
                elif nucleotide_2 == "G":
                    nuc_freq_2 = prob_g
                elif nucleotide_2 == "C":
                    nuc_freq_2 = prob_c
                
                for nucleotide_3 in self.nucleotides:
                    #codon = nucleotide_1 + nucleotide_2 + nucleotide_3
                    
                    if nucleotide_3 == "A":
                        nuc_freq_3 = prob_a
                    elif nucleotide_3 == "T":
                        nuc_freq_3 = prob_t
                    elif nucleotide_3 == "G":
                        nuc_freq_3 = prob_g
                    elif nucleotide_3 == "C":
                        nuc_freq_3 = prob_c
                    
                    frequency = nuc_freq_1 * nuc_freq_2 * nuc_freq_3
                    self.codon_dist[nucleotide_1][nucleotide_2][nucleotide_3] = round(frequency, 6)
                    
        return self.codon_dist

    def get_amino_acid_dist(self):
        """
        Return the expected amino acid distribution (fractions) based on the
        expected distribution (fractions) of the different codons.

        :return: Dictionary that contains the expected amino acid distribution.
                 The keys are the 20 different amino acids, the values are the
                 corresponding frequencies (rounded to 6 digits).
        """
        
        bases = "TCAG"
        #codons = [a + b + c for a in bases for b in bases for c in bases]
        self.amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
        self.cod_dist = [self.codon_dist[a][b][c] for a in bases for b in bases for c in bases] # computes distributions for codons
        self.amino_cod= list(zip(self.amino_acids, self.cod_dist)) # [ (F , 0.12334), (F , 0.23455)] aa & codon distr.
        self.aminotab = {}
        elem1 = 0 
        for elem in self.amino_cod:
            value = elem[1]    
            if elem[0] in self.aminotab: # if the amino acid key is already in the dict
                elem1 = self.aminotab[elem[0]] # store the old value
                self.aminotab[elem[0]] = elem1 + value # add new value to old vaue at the current key
                
            else:
                self.aminotab[elem[0]] = value
        
        normalization = sum(self.cod_dist) - (self.aminotab['*'])
        
        for k, v in self.aminotab.items():
            self.aminotab[k] = round((v / normalization), 6)
        
        del self.aminotab['*']
        
        
        
        
        
        return self.aminotab
