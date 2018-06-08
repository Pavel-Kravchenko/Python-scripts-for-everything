from sys import argv
indata = open(argv[1], "r")
line = indata.readline()
while line != "":
    if line[0] == ">":
        line = line[1:].split(" ")
        print line[0]
    line = indata.readline()    

indata.close()