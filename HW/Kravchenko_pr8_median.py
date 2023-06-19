from sys import argv

indata = open(argv[1], "r")

d = indata.readline()

d = d.strip()
d = d.split()
d = [float(i) for i in d]
f = sorted(d)
print f
g = len(f)
if g%2 != 0: 
    print "Median", f[g/2]
if g%2 == 0:
    h = (f[(len(f) - 1)/ 2] + f[(len(f) + 1)/ 2]) /2   
    print "Median", h 

indata.close()

