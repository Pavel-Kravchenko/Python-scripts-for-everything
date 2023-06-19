
# coding: utf-8

# In[8]:

import random


M = {0: 0, 1: 1}


def fibonacci(n):
    if n in M:
        return M[n]
    M[n] = fibonacci(n - 1) + fibonacci(n - 2)
    return M[n]


seed, min_lim, max_lim, count = map(int, input().strip().split())
random.seed(seed)
for i in range(count):
    fib_id = random.randint(min_lim, max_lim)
    fib = fibonacci(fib_id)
    print (fib)

