
# coding: utf-8

# In[32]:


from Bio import SeqIO

seql = []
f = open('in.txt', 'r')
lines = list(f.readlines())
for i in range(len(lines)):
    if lines[i] == '\n':
        mark = i
ids = lines[:mark]
with open('db.txt', 'w') as f:
    f.write(''.join(lines[(mark+1):]))
for j in ids:
    for i in SeqIO.parse("db.txt", "fasta"):
        if i.id == j.strip():
            seql.append(i)
SeqIO.write(seql, "out.txt", "fasta")

