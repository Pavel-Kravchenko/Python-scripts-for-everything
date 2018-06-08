#--------------biblio & modules-------------------------------------
from sys import argv
from Bio import AlignIO
from Bio.Alphabet import generic_nucleotide
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import GC
from Bio.Alphabet import generic_dna, generic_protein
from Bio.Data import CodonTable
from Bio.Seq import MutableSeq
from Bio.Seq import UnknownSeq
from Bio.Seq import reverse_complement, transcribe, back_transcribe, translate
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.SeqUtils.CheckSum import seguid
from random import randint
from Bio import Entrez



#--------------extracting the file-----------------------------
in_file = open(argv[1], "r")
out_file = open(argv[2], "w")

#--------------input fasta sequences-----------------------------------------------------
#for seq_record in SeqIO.parse(in_file, "fasta"):
    #print(seq_record.id)
    #print(repr(seq_record.seq))
    #print(len(seq_record))

#for record in SeqIO.parse(in_file, "fasta"):
     #print("%s %i" % (record.id, len(record)))

#--------------input genbank sequences-----------------------------------------------------
#for seq_record in SeqIO.parse(in_file,"genbank"):
    #print(seq_record.id)
    #print(repr(seq_record.seq))
    #print(len(seq_record))
    #print(seq_record.seq)





#---------------working with genbank sequences-----------------------------  
#records = list(SeqIO.parse(in_file,"genbank"))
#print(records[0].seq)  
#print(records[0].seq.count("AAA"))
#print(records[0].seq[0])
#print(records[-1].seq)
#print(len(records[0].seq))
#print(GC(records[0].seq))
#print(records[0].seq[4:12])
#print(records[0].seq[::-1])
#print(records[0].seq + records[1].seq) #for join files you should be sure that the alphabets are equal.





#----------------taking one seq from a file--------------------------------
#my_seq = Seq(str(records[0].seq), IUPAC.unambiguous_dna)
#print(my_seq)
#print(str(my_seq))
#print(my_seq.alphabet)

#record = SeqIO.read(in_file, "fasta")  #if one of one
#record = next(SeqIO.parse(in_file, "fasta")) #if one of many
#print("%s %i" % (record.id, len(record)))


#------------------------Input - Multiple Records--------------------------
#records = list(SeqIO.parse(in_file, "fasta"))
#print(len(records))
#print(records[1].id)

#record_dict = SeqIO.to_dict(SeqIO.parse(in_file, "fasta"))
#print(len(record_dict))
#print(record_dict["AY130325.1/1-585"])


#needs a filename directly (not a argv[])
#record_dict = SeqIO.index("in.fasta", "fasta")
#print(len(record_dict))
#print(record_dict.get_raw("AY130325.1/1-585").decode())
#print(record_dict["AY130325.1/1-585"].format("fasta"))
#record_dict.close()
#there is also fastQ format with the line end simvols





#----------------------Input - Alignments-------------------------------
#for record in SeqIO.parse(in_file, "fasta"):  #may be "clustal",..
     #print("%s %i" % (record.id, len(record)))
     #print(record.seq())





#---------------------output of records---------------------------------
#records = list(SeqIO.parse(in_file, "fasta"))
#SeqIO.write(records, out_file, "fasta")




#----------------placing seq in the table form-------------------------------------
#for index, letter in enumerate(my_seq):
    #print("%i %s" % (index, letter))






#------------script that gives a fasta format--------------------------
#fasta_format_string = ">" + records[0].id + "\n%s\n" % records[0].seq
#print(fasta_format_string)





#----------script that join sequences------------------------------
#list_of_seqs = [Seq(str(records[0].seq), generic_dna), Seq(str(records[1].seq), generic_dna), Seq(str(records[2].seq), generic_dna)]
#concatenated = Seq("", generic_dna)
#for s in list_of_seqs:
    #concatenated += s
#print(concatenated)
#You can use also ------sum(list_of_seqs, Seq("", generic_dna))-----------


#------------------format sequences--------------------------------- 
#print(records[0].seq.upper())   
#print(records[0].seq.lower())
#print(records[0].seq.complement())
#print(records[0].seq.reverse_complement())




#--------------------working with DNA RNA sequences---------------------------
#coding_dna = Seq("ATGGCCATTGTAAAAATCTGAUAAUAGUGA")
#back_transcript = messenger_rna.back_transcribe()
#full = template_dna.reverse_complement().transcribe()
#peptide_RNA = messenger_rna.translate()
#peptide_DNA = coding_dna.translate()
#Vertebrate_Mitochondrial_Code = coding_dna.translate(table=2)
#Yeast_Mitochondrial_Code  = coding_dna.translate(table=3)
#Yeast_Mitochondrial_Code_plus_stops  = coding_dna.translate(table=3, to_stop=True)
#Yeast_Mitochondrial_Code_plus_mail  = coding_dna.translate(table=3, stop_symbol="@")

#print(coding_dna)
#print(template_dna)
#print(messenger_rna)
#print(full)
#print(back_transcript)
#print(peptide_RNA)
#print(peptide_DNA)
#print(Vertebrate_Mitochondrial_Code)
#print(Yeast_Mitochondrial_Code)
#print(Yeast_Mitochondrial_Code_plus_stops)
#print(Yeast_Mitochondrial_Code_plus_mail)






#-------------------working with the full seq-----------------------------
#gene = Seq("GTGAAAAAGATGCAATCTATCGTACTCGCACTTTCCCTGGTTCTGGTCGCTCCCATGGCA" + \
#"AAGAAAGCTCCTCATGATCATCACGGCGGTCATGGTCCAGGCAAACATCACCGCTAA",generic_dna)
#protein = gene.translate()
#protein_bac = gene.translate(table="Bacterial")
#protein_to_stop = gene.translate(table="Bacterial", to_stop=True)
#protein_bac_cds = gene.translate(table="Bacterial", cds=True)

#print(gene)
#print(protein)
#print(protein_bac)
#print(protein_to_stop)
#print(protein_bac_cds)
    


#------------------Translation Tables--------------------------------
#standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
#or you can use ------- standard_table = CodonTable.unambiguous_dna_by_id[1]
#mito_table = CodonTable.unambiguous_dna_by_name["Vertebrate Mitochondrial"]
#or you can use ------- mito_table = CodonTable.unambiguous_dna_by_id[2]

#print(standard_table)
#print(mito_table.stop_codons)
#print(mito_table.start_codons)
#print(mito_table.forward_table["ACG"])




#--------------------------Comparing Seq objects---------------------------
#seq1 = Seq("ACGT", IUPAC.unambiguous_dna)
#seq2 = Seq("ACGT", IUPAC.ambiguous_dna)
#print(str(seq1) == str(seq2))

#dna_seq = Seq("ACGT", generic_dna)
#prot_seq = Seq("ACGT", generic_protein)
#print(dna_seq == prot_seq)
# it will giva the BiopythonWarning: Incompatible alphabets DNAAlphabet() and ProteinAlphabet() because of incorrect alph formats.




#---------------------creating MutableSeq-----------------
#my_seq = Seq("GCCATTGTAATGGGCCGCTGAAAGGGTGCCCGA", IUPAC.unambiguous_dna)
#mutable_seq = my_seq.tomutable()
#print(mutable_seq)

#mutable_seq = MutableSeq("GCCATTGTAATGGGCCGCTGAAAGGGTGCCCGA", IUPAC.unambiguous_dna)
#mutable_seq[5] = "C"
#mutable_seq.remove("T")
#mutable_seq.reverse()
#print(mutable_seq)
#you can edit or delete nucleotides as you want but there in no dict in this MutableSeq


#----------------------------UnknownSeq objects---------------------
#unk = UnknownSeq(20)
#unk_dna = UnknownSeq(20, alphabet=IUPAC.ambiguous_dna)
#print(unk)
#rint(unk_dna)


#--------------------------Working with strings directly---------------------
#my_string = "GCTGTTATGGGTCGTTGGAAGGGTGGTCGTGCTGCTGGTTAG"
#print(reverse_complement(my_string))
#print(back_transcribe(my_string))
#print(translate(my_string))


#Annotation!!!!!!!!!!!!!!!!!!!!!!!!!!Chapter!!!!!!!!!!!!!!!!!!!!






#An Example
#short_sequences = [] # Setup an empty list
#for record in SeqIO.parse(in_file, "fasta"):
#    if len(record.seq) < 590 :
       # Add this record to our list
#        short_sequences.append(record)

#print("Found %i sequences" % len(short_sequences))

#SeqIO.write(short_sequences, out_file, "fasta")








#-------------------futter--------------------------
in_file.close()
out_file.close()
print("")
print("#------------------------------------------------------------------#")
print(" -----------------------!!!Smile!!!--------------------------------")
print("                          The End")
