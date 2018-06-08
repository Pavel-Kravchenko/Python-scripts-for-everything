
# coding: utf-8

# In[ ]:

in_num = int(input())


def function(num):
    if num < 0:
        return 0
    if num < 40:
        return num
    out = function(num - 40) + function(num - 20)
    out = out + function(num - 5) + function(num - 1)
    return out


if in_num <= 10000:
    print(function(in_num))

