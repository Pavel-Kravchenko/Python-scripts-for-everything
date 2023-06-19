from sys import argv
import random
a = open("randomseq.fasta", "w")
b = int(argv[1])
a.write(">random")
w = ["A", "T", "G", "C"]
m = 0
l = ""
while m < b:
    t = str(random.choice(w)) 
    l = l + t
    m = m + 1  
X = 0
Y = 60
if len(l) > 60: 
    d = len(l)
    while d > 0:
        j = l[X:Y]
        d = d - 60
        a.write("\n" + j)
        X = X + 60
        Y = Y + 60      
else:
    a.write("\n" + l)