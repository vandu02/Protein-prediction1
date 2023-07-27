# -*- coding: utf-8 -*-
from Bio import SeqIO  # Tip: This module might be useful for parsing...


############ Exercise 3: SwissProt ##########
class SwissProt_Parser:
    PARSER = SeqIO

    def __init__(self, path, frmt='uniprot-xml'):
        """
            Initialize every SwissProt_Parser with a path to a XML-formatted UniProt file.
            An example file is included in the repository (P09616.xml).
            Tip: Store the parsed XML entry in an object variable instead of parsing it
            again & again ...
        """
        self.sp_parse = SeqIO.parse(path, frmt)
        self.sp_parse2 = SeqIO.parse(path,frmt)
        self.sp_anno = next(self.sp_parse) # Parse the XML file once and re-use it in the functions below

    # 2.2 SwissProt Identification
    def get_sp_identification(self):
        """
            Input:
                self: Use XML entry which has been parsed & saved during object initialization
            Return:
                identification: tuple consisting of in this order
                                    1. string: Unique SwissProt identifier for the given xml file
                                    2. string: Primary protein name
                                    3. string: Primary gene name
        """
        record = self.sp_anno
        id = record.id
        name = record.annotations["recommendedName_fullName"]
        gene = record.annotations["gene_name_primary"]
        identifier = (id, *name, gene)
        return identifier

    # 2.3 SwissProt Sequence Information
    def get_sp_sequence_info(self):
        """
            Input:
                self: Use XML entry which has been parsed & saved during object initialization
            Return:
                information: tuple consisting of in this order
                                    1. str: sequence of the UniProt entry
                                    2. int: sequence length of the UniProt entry
                                    3. int: sequence mass of the UniProt entry
        """
        seq_info = self.sp_anno
        prot = str(seq_info.seq)
        length = seq_info.annotations["sequence_length"]
        mass = seq_info.annotations["sequence_mass"]
        seq_len = (prot, length, mass)
        return seq_len

    # 2.4 Organism 
    def get_organism(self):
        """
            Input:
                self: Use XML entry which has been parsed & saved during object initialization
            Return:
                Return the name of the organsim as stated in the corresponding field
                of the XML data. Return value has to be a string.
        """
        seq_info = self.sp_anno
        organism = str(seq_info.annotations["organism"])
        return organism

    # 2.5 Localizations
    def get_localization(self):
        """
            Input:
                self: Use XML entry which has been parsed & saved during object initialization
            Return:
                Return the name of the subcellular localization as stated in the
                corresponding field.
                Return value has to be a list of strings.
        """
        seq_info = self.sp_anno
        localization = seq_info.annotations["comment_subcellularlocation_location"]
        return localization

    # 2.6 Cross-references to PDB
    def get_pdb_support(self):
        """
            Input:
                self: Use XML entry which has been parsed & saved during object initialization
            Return:
                Returns a list of all PDB IDs which support the annotation of the
                given SwissProt XML file. Return the PDB IDs as list.
        """
        pdb_ids = set()

        for record in self.sp_parse2:
            #print("hola")
            for db_reference in record.dbxrefs:
                if db_reference.startswith("PDB:"):
                    pdb_id = db_reference.split(":")[1]
                    pdb_ids.add(pdb_id)
                    
        print(pdb_ids)
        return list(pdb_ids)


def main():
    print('SwissProt XML Parser class')
    return None


if __name__ == '__main__':
    main()
