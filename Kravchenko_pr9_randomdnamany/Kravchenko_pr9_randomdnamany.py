from sys import argv
import random
fila = open("randomseqs.fasta", "w")
a = int(argv[1])
b = int(argv[2])
c = int(argv[3])
import mymodule

def tt(num):
    line = ss[num]
    return line

ss = []
if a > 0:
    for e in range(a):                  
        u = random.randint(b, c)        
        s = mymodule.i(u) 
        ss = ss + [s]                          

    for e in range(a):
        point = mymodule.ii(e)
        fila.write(point)
