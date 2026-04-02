
# coding: utf-8

# In[7]:


from Bio import SeqIO

sequences_list = []
sequences = SeqIO.parse("in.txt", "fasta")
for i in sequences:
    # print(i)
    i.seq = i.seq.complement()
    # print(i.seq)
    sequences_list.append(i)
SeqIO.write(sequences_list, "out.txt", "fasta")

