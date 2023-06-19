
# coding: utf-8

# In[ ]:

input_ = int(input())


def check2rec(num):
    if num == 1:
        return "YES"
    if num & 1:
        return "NO"
    return check2rec(num >> 1)


out_ = check2rec(input_)
print(out_)

