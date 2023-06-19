
# coding: utf-8

# In[41]:


from Bio import SeqIO

seql = []
seq = SeqIO.parse("in.txt", "fasta")
for i in seq:
    i.seq = i.seq.complement()
    seql.append(i)
SeqIO.write(seql, "out.txt", "fasta")

