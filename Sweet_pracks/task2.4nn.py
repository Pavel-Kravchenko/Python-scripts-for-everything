
# coding: utf-8

# In[11]:


input_n = int(input().strip())
_set_ = set()
for i in range(1, input_n + 1):
    _set_ = _set_.union(set(str(i)))
input_line = []

if input_n < 1000:
    while "HELP" not in input_line:
        input_line = input().split(' ')
        if "HELP" in input_line:
            break
        input_string_beatrice_2 = input()
        if input_string_beatrice_2 == "YES":
            YES = set(input_line)
            _set_ = _set_.intersection(YES)
        elif input_string_beatrice_2 == "NO":
            NO = set(input_line)
            _set_ = _set_.difference(NO)
_out_ = list(map(int, list(_set_)))
_out_ = sorted(fin)
for a in _out_:
    print(a, end=' ')

