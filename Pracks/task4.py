
# coding: utf-8

# In[16]:

line = input()
if len(line) <= 200:
    a = line.split("WUB")
    b = ""
    for i in a:
        if len(i) != 0:
            b = b + " " + i
print(b.strip())

