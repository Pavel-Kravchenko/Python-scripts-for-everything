from sys import argv
indata = open(argv[1], "r")
line = indata.readline()
line = indata.readline()
line = line.strip()
kol= 0
sux = 0
suy = 0
suz = 0
while len(line) != 0:
    lista = line.split("\t")
    sux = sux + float(lista[3])
    suy = suy + float(lista[4])
    suz = suz + float(lista[5])               
    line = indata.readline()         
    line = line.strip()              
    kol = kol + 1                    
print "The geometric center of the molacule is"
print "%.3f" % (sux / kol), "%.3f" % (suy / kol), "%.3f" % (suz / kol)          
    
