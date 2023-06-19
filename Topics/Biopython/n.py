
# coding: utf-8

# In[2]:


from Bio import SeqIO

seql = []
file = open('in.txt', 'r')
lines = list(file.readlines())
for i in range(len(lines)):
    if lines[i] == '\n':
        mark = i
list_id = lines[:mark]
with open('db.txt', 'w') as f:
    f.write(''.join(lines[(mark+1):]))
for j in list_id:
    for i in SeqIO.parse("db.txt", "fasta"):
        if i.id == j.strip():
            seql.append(i)
SeqIO.write(seql, "out.txt", "fasta")

