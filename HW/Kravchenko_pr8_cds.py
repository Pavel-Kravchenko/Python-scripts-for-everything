from sys import argv
indata = open(argv[1], "r")
line = indata.readline()
a = 0
while len(line) > 0:
    sline = line.strip()
    if len(sline) > 0:
        if "CDS" in sline:
           a = a + 1
           line = indata.readline()
    line = indata.readline()
print "There is ", a, "CDS in the file"
indata.close()