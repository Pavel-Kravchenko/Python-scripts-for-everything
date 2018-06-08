
# coding: utf-8

# In[45]:

with open('in.txt', 'r') as f:
    in_lines = f.readlines()
dic1 = dict()
dic2 = dict()
for i in in_lines:
    line = i.split()
    if line[0] == 'id1':
        ind = True
    if line[0][:2] != 'id' and ind is True:
        if line[0] not in dic1.keys():
            dic1[line[0]] = []
        dic1[line[0]].append(line[1])
    if line[0] == 'id2':
        ind = False
    if line[0][:2] != 'id' and ind is False:
        if line[0] not in dic2.keys():
            dic2[line[0]] = []
        dic2[line[0]].append(line[1])

ans = []
for i in dic1.keys():
    for j in dic1[i]:
        if j in dic2.keys():
            ans.append(i + '|' + ' '.join(dic2[j]))

fin = []
for i in range(len(ans)):
    ans[i] = ans[i].split('|')
    ans[i][1] = ans[i][1].split()
    for j in ans[i][1]:
        fin.append(ans[i][0] + ' ' + j)
fin.sort()

with open('out.txt', 'w') as f:
    f.write('id1 -> id3' + '\n')
    for i in fin:
        f.write(i + '\n')

