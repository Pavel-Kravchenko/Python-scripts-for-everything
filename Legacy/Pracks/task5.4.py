
# coding: utf-8

# In[19]:

import pickle
import json
import csv
dict_d = {}
with open('in.txt') as csvfile:
    iris_dataset = csv.reader(csvfile)
    header = next(iris_dataset)
    for i in header:
        dict_d[i] = []
    for row in iris_dataset:
        for i in range(len(row)):
            dict_d[header[i]].append(int(row[i]))
with open('out.txt', 'w') as out_file:
    out_file.write(json.dumps(dict_d, sort_keys=True, indent=4))


# In[ ]:



