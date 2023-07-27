def get_orfs(genome, min_num_aa):
    """
    Find and return all ORFs within the provided genome sequence.

    :param genome: String with the genome sequence.
    :param min_num_aa: Minimum number of amino acids in the translated protein
                       sequence.

    :return: List of tuples representing the ORFs (details in worksheet).
    """
    #1) Test to determine if genome is a valid DNA sequence
    """     nucleotides = ['A', 'C', 'G', 'T']
    for nucl in genome:
        if nucl not in nucleotides:
            raise TypeError("This is not a DNA") """
    verify(genome)
        
    #2) create  reading frames (first forward frames)
    
    reading_frames = [get_reading_frames(genome, i) for i in range (0,3)]
    reading_frames_rev = [get_reading_frames(reverse_compl(genome), i ) for i in range(0,3)]
    #print("reading frames", reading_frames)
    
    ORF_f1 = search_codon(reading_frames[0], 0, min_num_aa, False)
    ORF_f2 = search_codon(reading_frames[1], 1, min_num_aa, False)
    ORF_f3 = search_codon(reading_frames[2], 2, min_num_aa, False)
    fw_frames = [*ORF_f1, *ORF_f2, *ORF_f3] #list with all the forward frames
   
    ORF_rev1 = search_codon(reading_frames_rev[0], 0, min_num_aa, True)
    ORF_rev2 = search_codon(reading_frames_rev[1], 1, min_num_aa, True)
    ORF_rev3 = search_codon(reading_frames_rev[2], 2, min_num_aa, True)
    rev_frames = [*ORF_rev1, *ORF_rev2, *ORF_rev3] #list with all the reverse frames
    
    #2.1) find circular sequences
    
    #CIRC_fw = [find_circular_orfs(fw_frame, False) for fw_frame in fw_frames]
    #CIRC_rev = [find_circular_orfs(rev_frame, True) for rev_frame in rev_frames]
    
    length = len(genome) - 1
    
    CIRC_f0 = find_circular_orfs(reading_frames[0], False, length, min_num_aa)
    CIRC_f1 = find_circular_orfs(reading_frames[1], False, length, min_num_aa)
    CIRC_f2 = find_circular_orfs(reading_frames[2], False, length, min_num_aa)
    
    CIRC_r0 = find_circular_orfs(reading_frames_rev[0], True, length, min_num_aa)
    CIRC_r1 = find_circular_orfs(reading_frames_rev[1], True, length, min_num_aa)
    CIRC_r2 = find_circular_orfs(reading_frames_rev[2], True, length, min_num_aa)
    
    circles_fw = [CIRC_f0, CIRC_f1, CIRC_f2]
    #print(circles_fw)
    circles_rev = [CIRC_r0, CIRC_r1, CIRC_r2]
    
    #largest_circle_fw = max(circles_fw, key = circles_fw[2]) 
    largest_circle_fw = max(circles_fw, key=lambda x: x[2])
    largest_circle_rev = max(circles_rev, key=lambda x: x[2])
    #print(largest_circle_fw)
    #print(largest_circle_rev)
    
    #print("CIRC_F0: ", CIRC_f0, "\n")
    #print("CIRC_F1: ", CIRC_f1, "\n")
    #print("CIRC_F2: ", CIRC_f2, "\n")
    #print("CIRC_R0: ", CIRC_r0, "\n")
    #print("CIRC_R1: ", CIRC_r1, "\n")
    #print("CIRC_R2: ", CIRC_r2, "\n")
    #print("Circles forward: ", CIRC_fw)
    #print("Circles rev: ", CIRC_rev)
    
    #ORFseq_list.append([ORFsequence, init_pos, len(ORFsequence), aasequence, frame, currentCodon]) 
   
    #3) stop-codon grouping
    TAA = []
    TAG = []
    TGA = []
    
    for orf in fw_frames:
        stop = orf[5]
        if stop == 'TAA':
            TAA.append(orf)
        elif stop == 'TAG':
            TAG.append(orf)
        else:
            TGA.append(orf)
    for orf in rev_frames:
        stop = orf[5]
        if stop == 'TAA':
            TAA.append(orf)
        elif stop == 'TAG':
            TAG.append(orf)
        else:
            TGA.append(orf)
            
    #4) check for overlapping
    

    overlaps = is_overlapping(TAA, length)
    #print(overlaps)
    overlaps += is_overlapping(TGA, length)
    overlaps += is_overlapping(TAG, length)
    
    #5) Update positions for reverse frames:
    overlaps = get_reverse_pos(overlaps, length)


    #list_circular = [*largest_circle_fw]
    #list_circular.append([*largest_circle_rev])
    #res_fw = []
    #res_fw.append(largest_circle_fw)

    #res_rev = []
    #res_rev.append(largest_circle_rev)
    #overlaps += res_fw
    #overlaps += res_rev
    
    #6) Format output: [pos, length, protein, True/False]
    #circles = []
    #circles.append(largest_circle_fw)
    #circles.append(largest_circle_rev)
    #print("overlaps: ", overlaps)
    #print("circular: ", list_circular)
    
    overlaps2 = [(elem[1], elem[2], elem[3], elem[4]) for elem in overlaps]

    if largest_circle_fw[2] > 0:
      lcfw = (int(largest_circle_fw[1]), int(largest_circle_fw[2]), str(largest_circle_fw[3]), bool(largest_circle_fw[4]))
      overlaps2.insert(2, lcfw)

    if largest_circle_rev[2] > 0:
      lcrev = (int(largest_circle_rev[1]), int(largest_circle_rev[2]), str(largest_circle_rev[3]), bool(largest_circle_rev[4]))
      overlaps2.insert(3, lcrev)

    


    #output = overlaps2
    #output += circles_res
    #print('overlaps = ', overlaps)
    #return overlaps2
    return overlaps2
    
    #return unique_list
    #print(len(overlaps))





def reverse_compl(seq):
    ReverseCompDNA = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    rev_comp = ''.join([ReverseCompDNA[nuc] for nuc in seq])[::-1] #get reversed complement
    return rev_comp


def get_reading_frames (seq, init_pos):
    return ''.join([seq[pos] for pos in range(init_pos, len(seq), 1)])

# ORFseq_list.append([ORFsequence, init_pos, len(ORFsequence), aasequence, frame, currentCodon])
def is_overlapping(stop_group, length):
    results = []
    for i in stop_group:
        for j in stop_group[1:]:
            #if start_1 <= end_2 && start_2 <= end_1 (condition for overlapping)
            
            if i[4] == False and j[4] == False: #both are forward frames
                #print("Both are Forward")
                if i[1] <= (j[2] + j[1]) and j[1] <= (i[2] + i[1]): 
                    #print("OVERLAP")
                    #if they are overlapping, compare their lengths
                    if len(i[0]) > len(j[0]):
                        if j == stop_group[-1]:
                            results.append(i)
                            #print("i is the largest amongst all")
                else: #if not overlap
                    #print("no overlap")
                    if j == stop_group[-1]:
                        #save result i
                        results.append(i)
                        #print("i is the largest amongst all j")
                        
            elif i[4] == False and j[4] == True: #one is forward and the other reverse
                if i[1] <= length - (j[1] + j[2]) and length - j[1] <= (i[1] + i[2]): 
                    #print("OVERLAP")
                    #if they are overlapping, compare their lengths
                    if len(i[0]) > len(j[0]):
                        #print("i is larger than j")
                        if j == stop_group[-1]:
                            results.append(i)
                            #print("i is the largest amongst all")
                else: #if not overlap
                    #print("no overlap")
                    if j == stop_group[-1]:
                        #save result i
                        results.append(i)
                        #print("i is the largest amongst all j")
                        
            elif i[4] == True and j[4] == False: #one is reverse the other forward
                if length - i[1] <= (j[1] + j[2]) and j[1] <= length - (i[1] +i[2]): 
                    #print("OVERLAP")
                    #if they are overlapping, compare their lengths
                    if len(i[0]) > len(j[0]):
                        #print("i is larger than j")
                        if j == stop_group[-1]:
                            results.append(i)
                            #print("i is the largest amongst all")
                else: #if not overlap
                    #print("no overlap")
                    if j == stop_group[-1]:
                        #save result i
                        results.append(i)
                        #print("i is the largest amongst all j")
            
            elif i[4] == True and j[4] == True: #both are reverse
                if length - i[1] <= length - (j[1] + j[2]) and length - j[1] <= length - (i[1] + i[2]): 
                    #print("OVERLAP")
                    #if they are overlapping, compare their lengths
                    if len(i[0]) > len(j[0]):
                        #print("i is larger than j")
                        if j == stop_group[-1]:
                            results.append(i)
                            #print("i is the largest amongst all")
                else: #if not overlap
                    #print("no overlap")
                    if j == stop_group[-1]:
                        #save result i
                        results.append(i)
                        #print("i is the largest amongst all j")
                
    return results

def get_reverse_pos(overlaps, length):
    for elem in overlaps:
        if elem[4] == True:
            new_pos = length - elem[1]
            elem[1] = new_pos
            
    return overlaps
    
def search_codon(genome, idx, min_num_aas, frame):
    
    bases = "TCAG"
    codons = [a + b + c for a in bases for b in bases for c in bases]
    amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG' 
    gen_table = dict(zip(codons, amino_acids))
    STARTcodon = "ATG"
    STOPcodon = ["TAA", "TAG", "TGA"]
    position = 0
    ORFsequence = ""
    aasequence = ""
    ORFseq_list = []
    aaseq_list = []
  
    
    while (position < len(genome) - 2):
        # currentCodon-> take 3 steps.
        currentCodon = genome[position:position+3]
        num_aa = 0

        # currentCodon is True
        # currentCodon->codon='ATG'
        if (currentCodon == STARTcodon): 
            init_pos = position + idx # record position where START codon is being found
            while not (currentCodon in STOPcodon) and (currentCodon in gen_table):
                num_aa = num_aa + 1
                ORFsequence = ORFsequence + currentCodon
                aasequence += gen_table[currentCodon]
                position = position + 3 # Position for codon to move forward to select next codon
                # Replace currentCodon position
                currentCodon = genome[position:position+3]
            # Concatenate ORFs with respect to currentCodon position
            if num_aa >= min_num_aas: 
                ORFsequence = ORFsequence + currentCodon #once a stop codon has been found, add it to the ORF
                ORFseq_list.append([ORFsequence, init_pos, len(ORFsequence), aasequence, frame, currentCodon])
                #ORF_init.append(init_pos)
                aaseq_list.append(aasequence)
            #empty lists to start new reading:
            ORFsequence = ""
            aasequence = ""
        # Increase position
        position = position + 3
    return (ORFseq_list)

def verify (genome):
    nucleotides = ['A', 'C', 'G', 'T']
    for nucl in genome:
        if nucl not in nucleotides:
            raise TypeError("This is not a DNA")
    return True

def find_circular_orfs(genome, frame, length, min_num_aa):
    start_codon = 'ATG'
    stop_codons = ['TAA', 'TAG', 'TGA']
    bases = "TCAG"
    codons = [a + b + c for a in bases for b in bases for c in bases]
    amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG' 
    gen_table = dict(zip(codons, amino_acids))
    
   
    flag1 = True
    flag2 = True
    flag3 = True
    position = len(genome) - 3  
    position_end = len(genome) - 3
    position_lstop = 0 #start at the beginning of the string

    #Find rightmost END Codon
    
    while (flag2 and position_end >= 0):
        # currentCodon-> take 3 steps.
        currentStop = genome[position_end:position_end + 3]
        #print("currentStop Loop", currentStop)
        
        # STOPCodon is True
       
        if (currentStop in stop_codons): 
            stop_codon_pos = position_end  # record position where START codon is being found
            flag2 = False
          #  print("Stop Codon pos", stop_codon_pos)
        else:
            position_end = position_end - 3
     
    
    #Find START codon at the right of the rightmost stop codon
    
    right_start_pos = stop_codon_pos
    
    while (flag1 and right_start_pos <= len(genome) - 3 ):
        
        rightStartCodon = genome[right_start_pos: right_start_pos + 3]
        if (rightStartCodon == start_codon):

            flag1 = False
        else:
            right_start_pos = right_start_pos + 3
            
        
    
    #Find leftmost STOP codon, traversing from left to right
    
    while (flag3 and position_lstop <= len(genome) - 3):
        # currentCodon-> take 3 steps.
        currentLStop = genome[position_lstop:position_lstop + 3]

        # STOPCodon is True
       
        if (currentLStop in stop_codons): 
            left_stop_codon_pos = position_lstop  # record position where START codon is being found
            flag3 = False
        else:
            position_lstop = position_lstop + 3

    #print("Rightmost STOP", currentStop, position_end)
    #print("position of START after rightmost stop = ",rightStartCodon, right_start_pos)
    #print("Leftmost STOP", currentLStop, left_stop_codon_pos, "\n")
   
    circular_orf = genome[right_start_pos: ] + genome[:left_stop_codon_pos + 3]
    
    #compute amino acid sequence
    protein = ""
    for i in range(0, len(circular_orf) - 5, 3 ):
        protein += gen_table[circular_orf[i:i+3]]
        
    # if circular dna is found in a reversed frame, then recalculate position:
    if frame == True:
        right_start_pos = length - right_start_pos
    
    #if amino acid sequence has at least min_num_aa:
    
    if len(protein) >= min_num_aa:
        return ([circular_orf, right_start_pos, len(circular_orf), protein, frame])
    else:
        return([0,0,0,0,frame])

    
    

