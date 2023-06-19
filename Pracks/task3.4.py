
# coding: utf-8

# In[7]:

in_num = int(input())
M = {0: 0, 1: 1}


def function(num):
    if num < 0:
        return 0
    if num < 40:
        return num
    if num in M:
        return M[num]
    out = function(num - 40) + function(num - 20)
    out = out + function(num - 5) + function(num - 1)
    M[num] = out
    return out


if in_num <= 100000:
    print(function(in_num))

