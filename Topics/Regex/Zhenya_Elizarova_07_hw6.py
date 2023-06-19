
# coding: utf-8

# In[5]:


import re


text = input()
l=[]
for match in re.finditer(r'[+-]?([0-9]*)?[\s|.][0-9]+\s', text):
    l.append(match.group(0))
for i in l:
    text1=re.sub(i, ' '+str(float(i)**3)+' ', text)
    text=text1
print(text1)

