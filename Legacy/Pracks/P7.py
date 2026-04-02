
# coding: utf-8

# In[180]:


f = open('in.txt')
file_content = f.readlines()
f.close()
a = dict()
b = dict()
ind = 0
for i in file_content:
    if i.split()[0] == 'id2':
        ind = 1
    if i.split()[0] != 'id1' and ind == 0:
        if i.split()[0] not in a.keys():
            a[i.split()[0]] = []
        a[i.split()[0]].append(i.split()[1])
    if i.split()[0] != 'id2' and ind == 1:
        if i.split()[0] not in b.keys():
            b[i.split()[0]] = []
        b[i.split()[0]].append(i.split()[1])
ans = []
for i in a.keys():
    for j in a[i]:
        if j in b.keys():
            ans.append(i + '.' + ' '.join(b.get(j)))
fin = []
for i in range(len(ans)):
    ans[i] = ans[i].split('.')
    ans[i][1] = ans[i][1].split()
    for j in ans[i][1]:
        fin.append(ans[i][0] + ' ' + j)
fin.sort()
with open('out.txt', 'w') as f:
    f.write('id1 -> id3')
    for i in fin:
        f.write('\n' + i)

