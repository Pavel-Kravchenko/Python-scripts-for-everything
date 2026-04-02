
# coding: utf-8

# In[48]:

import itertools
in_num = int(input().strip())
if 0 < in_num < 10:
    in_num = [i for i in range(0, in_num)]
    if len(in_num) == 1:
        print(in_num)
    else:
        out = list(itertools.permutations(in_num))
        out = sorted(out)
        for i in range(len(out)):
            if len(out[i]) != 1:
                print(list(out[i]))

