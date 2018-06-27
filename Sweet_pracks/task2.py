
# coding: utf-8

# In[129]:


from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

calc = 1
for rec in SeqIO.parse("in.txt", "genbank"):
    with open("out.txt", "w") as f:
        if rec.features:
            for feature in rec.features:
                strand = feature.location.strand
                if feature.type == "CDS":
                    coding_genes = (
                        SeqRecord(seq=feature.location.extract(rec).seq,
                                  id=str(calc), description=str(strand)))
                    SeqIO.write(coding_genes, f, "fasta")
                    calc += 1

