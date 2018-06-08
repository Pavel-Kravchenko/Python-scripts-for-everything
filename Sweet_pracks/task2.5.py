
# coding: utf-8

# In[41]:

import operator
_input_ = input()
array = list(_input_.strip().split())
result = {i: array.count(i) for i in array}
result = sorted(result.items(), key=operator.itemgetter(0))
for i in result:
    print(str(i[0]) + " " + str(i[1]))

