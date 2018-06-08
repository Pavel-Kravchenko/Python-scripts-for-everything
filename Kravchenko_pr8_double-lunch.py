from sys import argv
indata = open(argv[1], "r")
line = indata.readline()
a = []
y = []
print "Wanted"
print "They are very dishonest:"
sline = line.strip()
while sline != "STOP":
    a = a + [sline]
    line = indata.readline()
    sline = line.strip() 

if a == []:
    print "None"
for r in a: 
    if r not in y:
        y = y + [r] 
        d = a.count(r)    
        if d > 1:         
            print r       
indata.close()
               
