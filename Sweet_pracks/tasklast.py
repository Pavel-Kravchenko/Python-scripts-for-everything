
# coding: utf-8

# In[47]:

input_n = int(input().strip())


def gen(n: int, counter_open: int, counter_close: int, ans: ""):
    if counter_open + counter_close == 2 * n:
        if counter_open + counter_close != 0:
            print(ans)
            return
    if counter_open < n:
        gen(n, counter_open + 1, counter_close, ans + '(')
    if counter_open > counter_close:
        gen(n, counter_open, counter_close + 1, ans + ')')


a = gen(input_n, 0, 0, "")

