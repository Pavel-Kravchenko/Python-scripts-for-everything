
# coding: utf-8

# In[26]:

line = input()
ide = "hello"
c = 0
if len(line) >= 1 and len(line) <= 100:
    for i in line:
        if i == ide[c]:
            c += 1
            if c >= 5:
                print("YES")
                break   
if c < 5: 
    print("NO")

