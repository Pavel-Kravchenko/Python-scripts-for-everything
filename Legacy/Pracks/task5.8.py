
# coding: utf-8

# In[3]:

with open('in.txt', 'r') as f:
    in_lines = f.readlines()
dic1 = dict()
dic2 = dict()
ind = 0
for i in in_lines:
    line = i.split()
    if ind == 1 and len(line) != 3:
        dic1[line[0]] = float(line[1])
    if ind == 2 and len(line) != 3:
        dic2[line[0]] = float(line[1])
    if len(line) == 3:
        ind += 1
ans = []
for i in dic1.keys():
    if i in dic2.keys():
        ans.append(i + ' ' + str(round((dic1.get(i) + dic2.get(i)), 5)))
    if i not in dic2.keys():
        ans.append(i + ' NODATA')
for j in dic2.keys():
    if j not in dic1.keys():
        ans.append(j + ' NODATA')
with open('out.txt', 'w') as f:
    f.write('id -> sum' + '\n')
    for i in sorted(ans):
        f.write(i + '\n')

