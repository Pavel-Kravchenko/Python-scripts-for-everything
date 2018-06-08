
# coding: utf-8

# In[23]:

import random


def fib(n):
    v1, v2, v3 = 1, 1, 0
    for rec in bin(n)[3:]:
        calc = v2 * v2
        v1, v2, v3 = v1 * v1 + calc, (v1 + v3) * v2, calc + v3 * v3
        if rec == '1':
            v1, v2, v3 = v1 + v2, v1, v2
    return v2


seed, min_lim, max_lim, count = map(int, input().strip().split())
random.seed(seed)
for i in range(count):
    fib_id = random.randint(min_lim, max_lim)
    fib1 = fib(fib_id)
    print(fib1)

