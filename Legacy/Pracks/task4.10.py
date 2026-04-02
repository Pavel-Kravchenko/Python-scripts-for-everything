
# coding: utf-8

# In[112]:


from Bio import SeqIO

lines = []
with open('in.txt', 'r') as f:
    line = f.readline().strip()
    while len(line) != 0:
        lines.append(line)
        line = f.readline().strip()
records = SeqIO.parse("in.txt", "fasta")
sequence_list = []
for i in lines:
    for j in records:
        if j.id == i:
            sequence_list.append(j)
            break
        else:
            pass
SeqIO.write(sequence_list, "out.txt", "fasta")

