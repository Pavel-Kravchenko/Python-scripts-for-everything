
# coding: utf-8

# In[2]:

import functools
import sys
import threading
sys.setrecursionlimit(100000)
threading.stack_size(0x2000000)


def fibonacci_rec_hash(n, calc_memory={0: 0, 1: 1}):
    if n < 0:
        return 0
    result = calc_memory.get(n, None)  # safe way to get value,
    # if value doesn't exist, returns default value, here it's None
    if result is not None:
        return result

    result = fibonacci_rec_hash(n - 1) + fibonacci_rec_hash(n - 2)
    calc_memory[n] = result
    return result


in_num = int(input())
M = {0: 0, 1: 1}


@functools.lru_cache(maxsize=1000, typed=False)
def function(num):
    if num < 0:
        return 0
    if num < 40:
        return num
    if num in M:
        return M[num]
    out = function(num - 40) + function(num - 20)
    out = out + function(num - 5) + function(num - 1)
    M[num] = out
    return out


if in_num <= 100000:
    print(function(in_num))

