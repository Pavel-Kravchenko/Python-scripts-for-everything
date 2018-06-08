import math 
from sys import argv
a = float(argv[1])
 
a = math.radians(a)


z = math.sin(a)
x = math.cos(a)

                
print format(z, ".4f")
print format(x, ".4f")
        


