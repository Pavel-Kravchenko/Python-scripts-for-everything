
# coding: utf-8

# In[21]:

import pickle
import json
dic = pickle.load(open("in", "rb"))
with open('out.txt', 'w') as out_file:
    out_file.write(json.dumps(dic, sort_keys=True, indent=4))

