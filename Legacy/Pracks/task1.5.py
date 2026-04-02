
# coding: utf-8

# In[14]:

with open('in.txt', "r") as in_file:
    lenl_all = in_file.readlines()
    counter = 0
    in_a = 0
    for seq in lenl_all:
        lenl = seq.strip()
        a_in_line = list(lenl).count("a")
        in_a = in_a + a_in_line
        counter += 1
with open('out.txt', "w") as out_file:
    out_file.write(str(float(in_a/counter)))


# In[ ]:



