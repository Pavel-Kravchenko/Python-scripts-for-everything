
# coding: utf-8

# In[62]:

in_1 = input().strip().split(" ")
M = {0: 0, 1: 1}


def fib(n):
    if n in M:
        return M[n]
    M[n] = fib(n - 1) + fib(n - 2)
    return M[n]


bim_bam_bom = []
for i in range(0, int(in_1[1])):
    bim_bam_bom = bim_bam_bom + [input().strip()]
for i in range(len(bim_bam_bom)):
    print(bim_bam_bom[len(bim_bam_bom) - i - 1])
print(fib(int(in_1[0])))
for i in bim_bam_bom:
    print(i)

