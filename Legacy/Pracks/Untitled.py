
# coding: utf-8

# In[2]:


max_integer = int(input())

numbers = []
whole = set()
for m in range(1, max_integer+1):
    whole.add(str(m))

while 'HELP' not in numbers:
    numbers = input().split(' ')
    if 'HELP' in numbers:
        break
    question = input()
    if question == 'YES':
        yes = set(numbers)
        whole = whole.intersection(yes)

    elif question == 'NO':
        no = set(numbers)
        whole = whole.difference(no)

fin = list(map(int, list(whole)))
fin = sorted(fin)
for a in fin:
    print(a, end=' ')

