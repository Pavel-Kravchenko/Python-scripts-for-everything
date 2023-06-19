
# coding: utf-8

# In[28]:


import re


a=input()
num1=a.split()
regular = r"[АВЕКМНОРСТУХ]{1}\d{3}[АВЕКМНОРСТУХ]{2}\d{2,3}" 
regular1 = r"[АВЕКМНОРСТУХ]{2}\d{4,5}"
pat = re.search(regular, str(num1))
pat1 = re.search(regular1, str(num1))
if pat!= None:
    print('Private')
elif pat1!=None:
    print('Taxi')
else:
    print('Fail')


# In[108]:


import re

text = input() 
pattern=r'([0-1][0-9]|[2][1-3]):([0-5][0-9])(:[0-5][0-9])?'
pat = re.search(pattern, str(text))
print(re.sub(pattern, '(TBD)', text)) 


# In[96]:


import re


text = input()
for match in re.finditer(r'\d+(.\d+)?', text):
    l.append(match.group(0))
for i in l:
    text1=re.sub(i, str(float(i)**3), text)
    text=text1
print(text1)

