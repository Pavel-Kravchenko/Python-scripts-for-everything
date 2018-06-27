
# coding: utf-8

# In[36]:


from Bio import SeqIO

record_dict = SeqIO.to_dict(SeqIO.parse("in.txt", "fasta"))
#  print(record_dict)
new_dict = record_dict.copy()
for i in record_dict.keys():
    new_seq = record_dict[i].seq.complement()
    new_dict[i].seq = new_seq
#  print(new_dict)
SeqIO.write(new_dict.values(), "out.txt", "fasta")

