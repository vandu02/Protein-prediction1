# -*- coding: utf-8 -*-

from Bio.PDB.MMCIFParser import MMCIFParser  # Tip: This module might be useful for parsing...
import numpy as np
import Bio


############# Exercise 2: Protein Data Bank #############
# General remark: In our exercise every structure will have EXACTLY ONE model.
# This is true for nearly all X-Ray structures. NMR structures have several models.
class PDB_Parser:
    CIF_PARSER = MMCIFParser()  # parser object for reading in structure in CIF format

    def __init__(self, path):
        """
            Initialize every PDB_Parser with a path to a structure-file in CIF format.
            An example file is included in the repository (7ahl.cif).
            Tip: Store the parsed structure in an object variable instead of parsing it
            again & again ...
        """
        self.structure = self.CIF_PARSER.get_structure("struct", path)  # Parse the structure once and re-use it in the functions below

    # 2.8 Chains    
    def get_number_of_chains(self):
        """
            Input:
                self: Use Biopython.PDB structure which has been stored in an object variable
            Return:
                Number of chains in this structure as integer.
        """
        chains = list(self.structure.get_chains())
        n_chains = len(chains)
        return n_chains

    # 2.9 Sequence  
    def get_sequence(self, chain_id):
        """
            Input:
                self: Use Biopython.PDB structure which has been stored in an object variable
                chain_id  : String (usually in ['A','B', 'C' ...]. The number of chains
                        depends on the specific protein and the resulting structure)
            Return:
                Return the amino acid sequence (single-letter alphabet!) of a given chain (chain_id)
                in a Biopython.PDB structure as a string.
        """
        chain = self.structure[0][chain_id] #get chains
        residues = chain.get_residues() # get residues from chains

        sequence = ""
        for residue in residues:
            if residue.get_resname() in Bio.PDB.Polypeptide.standard_aa_names:
                sequence += Bio.PDB.Polypeptide.three_to_one(residue.get_resname())

        #sequence = 'SEQWENCE'
        return sequence

    # 2.10 Water molecules
    def get_number_of_water_molecules(self, chain_id):
        """
            Input:
                self: Use Biopython.PDB structure which has been stored in an object variable
                chain_id  : String (usually in ['A','B', 'C' ...]. The number of chains
                        depends on the specific protein and the resulting structure)
            Return:
                Return the number of water molecules of a given chain (chain_id)
                in a Biopython.PDB structure as an integer.
        """
        n_waters = 0
        for mod in self.structure:
            for chain in mod:
                if chain.id == chain_id:
                    for residue in chain:
                        if residue.get_resname().strip() == "HOH":
                            n_waters += 1


        #n_waters = 12
        return n_waters

    # 2.11 C-alpha distance
    def get_ca_distance(self, chain_id_1, index_1, chain_id_2, index_2):
        """
            Input:
                self: Use Biopython.PDB structure which has been stored in an object variable
                chain_id_1 : String (usually in ['A','B', 'C' ...]. The number of chains
                                depends on the specific protein and the resulting structure)
                index_1    : index of a residue in a given chain in a Biopython.PDB structure
                chain_id_2 : String (usually in ['A','B', 'C' ...]. The number of chains
                            depends on the specific protein and the resulting structure)
                index_2    : index of a residue in a given chain in a Biopython.PDB structure

                chain_id_1 and index_1 describe precisely one residue in a PDB structure,
                chain_id_2 and index_2 describe the second residue.

            Return:
                Return the C-alpha (!) distance between the two residues, described by
                chain_id_1/index_1 and chain_id_2/index_2. Round the returned value via int().

            The reason for using two different chains as an input is that also the distance
            between residues of different chains can be interesting.
            Different chains in a PDB structure can either occur between two different proteins
            (Heterodimers) or between different copies of the same protein (Homodimers).
        """
        residue_1 = self.structure[0][chain_id_1][index_1]
        residue_2 = self.structure[0][chain_id_2][index_2]
        atom_1 = residue_1["CA"]        
        atom_2 = residue_2["CA"]

        ca_distance = atom_1 - atom_2

        #print(int(ca_distance))



        #ca_distance = 12.56
        return int(ca_distance)

    # 2.12 Contact Map  
    def get_contact_map(self, chain_id):
        """
            Input:
                self: Use Biopython.PDB structure which has been stored in an object variable
                chain_id  : String (usually in ['A','B', 'C' ...]. The number of chains
                        depends on the specific protein and the resulting structure)
            Return:
                Return a complete contact map (see description in exercise sheet)
                for a given chain in a Biopython.PDB structure as numpy array.
                The values in the matrix describe the c-alpha distance between all residues
                in a chain of a Biopython.PDB structure.
                Only integer values of the distance have to be given (see below).
        """

        #length = 10

        chain = self.structure[0][chain_id]
        residues = [residue for residue in chain.get_residues() if "CA" in residue]
        length = len(list(residues))

        contact_map = np.zeros((length, length), dtype=np.float32)
                    
        for i, residue_i in enumerate(residues):
                for j, residue_j in enumerate(residues):
                    atom_i = residue_i["CA"].get_coord()
                    atom_j = residue_j["CA"].get_coord()
                    ca_distance = np.linalg.norm(atom_i - atom_j)
                    contact_map[i, j] = int(ca_distance)

        
        return contact_map.astype(np.int64)  # return rounded (integer) values
        

        """   
            for i in range(length):
            for j in range(i+1, length):
                distance = self.get_ca_distance(chain_id, i, chain_id, j)
                contact_map[i, j] = int(distance)
                contact_map[j, i] = int(distance)
        """

        #return contact_map.astype(np.int64)  # return rounded (integer) values

    # 2.13 B-Factors    
    def get_bfactors(self, chain_id):
        """
            Input:
                self: Use Biopython.PDB structure which has been stored in an object variable
                chain_id  : String (usually in ['A','B', 'C' ...]. The number of chains
                        depends on the specific protein and the resulting structure)
            Return:
                Return the B-Factors for all residues in a chain of a Biopython.PDB structure.
                The B-Factors describe the mobility of an atom or a residue.
                In a Biopython.PDB structure B-Factors are given for each atom in a residue.
                Calculate the mean B-Factor for a residue by averaging over the B-Factor
                of all atoms in a residue.
                Sometimes B-Factors are not available for a certain residue;
                (e.g. the residue was not resolved); insert np.nan for those cases.

                Finally normalize your B-Factors using Standard scores (zero mean, unit variance).
                You have to use np.nanmean, np.nanvar etc. if you have nan values in your array.
                The returned data structure has to be a numpy array rounded again to integer.
        """
        chain = self.structure[0][chain_id]
        residues = [residue for residue in chain.get_residues() if "CA" in residue]
        length = len(residues)
        
        #b_factors = np.array(length, dtype=np.float32)
        b_factors = np.empty(length, dtype=np.float32)


        for i, residue in enumerate(residues):
            atom_b_factors = [atom.get_bfactor() for atom in residue if atom.get_bfactor() is not None]
            if len(atom_b_factors) > 0:
                residue_b_factor = np.nanmean(atom_b_factors)
            else:
                residue_b_factor = np.nan

            b_factors[i] = residue_b_factor

        # Normalize B-Factors using standard scores
        mean = np.nanmean(b_factors)
        std = np.nanstd(b_factors)
        b_factors = (b_factors - mean) / std
        
       
        return b_factors.astype(np.int64)  # return rounded (integer) values


def main():
    print('PDB parser class.')
    return None


if __name__ == '__main__':
    main()
