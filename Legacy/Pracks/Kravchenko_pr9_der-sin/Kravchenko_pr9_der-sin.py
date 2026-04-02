import math 
from sys import argv
x = float(argv[1])
d = float(argv[2]) 
f = math.sin(x)
g = (math.sin(x + d) - math.sin(x - d))/(2*d)

p = math.cos(x)
l = abs(p - g)                
print str(format(x, ".4f")) + "\t" + str(format(g, ".4f")) + "\t" + str(format(p, ".4f")) + "\t" + str(format(l, ".4f"))
