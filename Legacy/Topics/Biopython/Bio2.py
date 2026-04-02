
# coding: utf-8

# In[ ]:


from Bio import AlignIO

index = []
alignment = AlignIO.read('in.txt', "stockholm")
for i in range(len(alignment[0])):
    if '-' not in alignment[:, i]:
        index.append(i)
result = alignment[:, index[0]:(index[0] + 1)]
for j in index[1:]:
    result = result + alignment[:, j:(j+1)]
AlignIO.write(result, 'out.txt', "clustal")

