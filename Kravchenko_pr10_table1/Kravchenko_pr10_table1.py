from sys import argv
indata = open(argv[1], "r")
line = indata.readline()
line = indata.readline()
line = line.strip()
kol= 1
pu = []
while len(line) != 0:
    l = line[4]
    if l not in ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "V", "U", "W", "X", "Y", "Z"]:
        print "invalid data" 
        raise SystemExit(1)
    if l not in pu:          
        pu = pu + [l]                
    line = indata.readline()         
    line = line.strip()              
    kol = kol + 1                    
            
print kol - 1, "residues"
print len(pu), "chain"
    
