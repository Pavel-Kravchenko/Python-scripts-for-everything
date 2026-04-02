
# coding: utf-8

# In[14]:


def lcm(a, b):
    m = a * b
    while a != 0 and b != 0:
        if a > b:
            a %= b
        else:
            b %= a
    return m // (a + b)


input_l = input().strip().split()
one = int(input_l[0])
two = int(input_l[1])
print(lcm(one, two))

