
# coding: utf-8

# In[2]:


n_lines = int(input())
names = {}
if n_lines <= 10000 and n_lines >= 1:
    for i in range(0, n_lines):
        user = input().strip()
        if user not in names:
            names[user] = 1
            print('OK')
        else:
            print(user + str(names.get(user)))
            names[user] += 1

