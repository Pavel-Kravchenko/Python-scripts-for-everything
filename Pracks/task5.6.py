
# coding: utf-8

# In[19]:

with open('in.pdb', "r") as in_file:
    lenl_all = in_file.readlines()
    for seq in lenl_all:
        seq = seq.split()
        if len(seq) >= 3 and seq[2] == "HOH" and seq[0] == "FORMUL":
            stri = seq[3][1:]
            strin = ""
            for i in stri:
                if i != "(":
                    strin = strin + i
                else:
                    break
with open('out.txt', "w") as out_file:
    out_file.write(strin)

