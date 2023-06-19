
# coding: utf-8

# In[44]:

a = input()
a = a.split()
c = 1
b = 0
for i in a:
    if int(i) <= 109 and int(i) >= 1:
        if i in a[c:]:
            b += 1
        c += 1
print(b)

