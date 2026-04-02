from sys import argv
indata = open(argv[1], "r")
c = argv[2]
line = indata.readline()
line = line.strip()    
while len(line) > 0:
    line = line.strip()
    if c in line:
        print line
    line = indata.readline()
indata.close()
 