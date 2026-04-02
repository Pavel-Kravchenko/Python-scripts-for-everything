
# coding: utf-8

# In[198]:


f = open('in.txt')
file_content = f.readlines()
f.close()
a = dict()
b = dict()
ind = 0
for i in file_content:
    if ind == 1 and len(i.split()) != 3:
        a[i.split()[0]] = float(i.split()[1])
    if ind == 2 and len(i.split()) != 3:
        b[i.split()[0]] = float(i.split()[1])
    if len(i.split()) == 3:
        ind += 1
ans = []
for i in a.keys():
    if i in b.keys():
        ans.append(i + ' ' + str(round((a.get(i) + b.get(i)), 5)))
    if i not in b.keys():
        ans.append(i + ' NODATA')
for j in b.keys():
    if j not in a.keys():
        ans.append(j + ' NODATA')
with open('out.txt', 'w') as f:
    f.write('id -> sum')
    for i in sorted(ans):
        f.write('\n' + i)

