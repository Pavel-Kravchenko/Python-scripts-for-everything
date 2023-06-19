
# coding: utf-8

# In[19]:

from decimal import *
getcontext().prec = 100
with open('in.txt', "r") as in_file:
    lenl_all = in_file.readlines()
    summ = 0
    for num in lenl_all:
        if len(num.strip()) != 0:
            try:
                num = float(num.strip())
                summ = Decimal(summ) + Decimal(num)
            except ValueError:
                summ = "nan"
                break
with open('out.txt', "w") as out_file:
    out_file.write(str(summ))

