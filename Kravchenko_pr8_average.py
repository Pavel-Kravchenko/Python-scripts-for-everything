from sys import argv

indata = open(argv[1], "r")

d = indata.readline()

d = d.strip()
d = d.split()
d = [float(i) for i in d]
o = 0
s = 0
while len(d[o:]) > 0:
    
    s = s + d[o]
    o = o + 1 
 
print "Average is", s/o

indata.close()

