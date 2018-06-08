from sys import argv
import urllib2
a = urllib2.urlopen("https://files.rcsb.org/download/" + str(argv[1]) + ".pdb")
line = a.readline()
name = argv[1] + ".tab"
fila = open(name, "w")
sd = ["Residue", "Chain", "Number", "X", "Y", "Z"]
sd = "\t".join(sd)
fila.write(sd + "\n")
import reader
while len(line) > 0:
    retrouver = reader.wolf(line)
    if len(retrouver) != 0:
        fila.write(retrouver + "\n")
    line = a.readline()    
    

 