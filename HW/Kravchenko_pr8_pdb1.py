from sys import argv
indata = open(argv[1], "r")
c = argv[2]
line = indata.readline()
n = 0    
amc = c.upper()
while len(line) > 0:
    line = line.strip()
    if len(line) > 0:
        if "ATOM  " in line:
            if (amc) in line[17:20]:
                if " CA " in line[12:16]:   
                    n = n + 1   
    line = indata.readline()
print "The number of your amino acides is "
print n
indata.close()
 