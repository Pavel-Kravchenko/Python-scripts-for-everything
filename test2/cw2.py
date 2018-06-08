from sys import argv
indata = open(argv[1], "r")
line = indata.readline()
t = []
import myfunc
while len(line) != 0:
    l = myfunc.func(line)  
    if l != "":
        t = t + [l]
    line = indata.readline()
for p in t:
    print p