from sys import argv 
c = int(argv[1])
import fb
sas = []
for i in range(c + 1):
    a = fb.fib(i) 
    sas =  sas + [str(a)] 
     
print " ".join(sas)
 
