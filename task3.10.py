
# coding: utf-8

# In[28]:


from Bio import AlignIO
from Bio import SeqIO


def cut_a_gap(align):
    n = float(len(align[0]))
    i = 0
    while i < n:
        if align[:, i].count('-') != 0:
            if i == 0:
                align = align[:, 1:]
            elif i+1 == n:
                align = align[:, :i]
            else:
                align = align[:, :i] + align[:, i+1:]
            n -= 1  # seq. 1 shorter
        else:    # nothing to delete, proceed
            i += 1
    return align


def cut_gap_in_blocks():   # usefull for the future, not for this script
    n = float(len(align[0]))
    i = 0
    while i < n:
        ct = 0
        while i+ct < n and align[:, i+ct].count('-') / n > 0.5:
            ct += 1

        if ct > 0:       # delete columns [i:i+ct]
            if i == 0:
                align = align[:, ct:]
            elif i+ct == n:
                align = align[:, :i]
            else:
                align = align[:, :i] + align[:, i+ct:]
            n -= ct    # seq. ct positions shorter
        else:    # nothing to delete, proceed
            i += 1
        return align


alignment = AlignIO.read("in.txt", "stockholm")
# print(alignment)
out = cut_a_gap(alignment)

SeqIO.write(out, "out.txt", "clustal")

