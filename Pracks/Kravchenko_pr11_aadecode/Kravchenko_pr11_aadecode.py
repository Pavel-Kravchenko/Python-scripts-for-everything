from sys import argv
import urllib2
a = urllib2.urlopen("http://kodomo.fbb.msu.ru/FBB/BasicData/aa_properties/aacodes.txt")
line = a.readline()
line = a.readline()
line = line.strip()
stri = line.split()
dica = {}
dicaaa = {}
                               
while len(line) > 0:
    dica[stri[0]] = stri[1]  
    dicaaa[stri[1]] = stri[0]                       
    line = a.readline()        
    line = line.strip()        
    stri = line.split()

def func(dic, data):
    if data in dic.keys():
        return dic[data]
    if data not in dic.keys():
        return 0 
if len(argv) == 1:
    a = ""
    argv = []
    while len(a) != 1 and len(a) != 3:
        a = str(raw_input("Something is wrong. The program can not use your data. Please try to write them in other way."))
    argv = argv + [""] +[a]

if len(argv[1]) != 1 and len(argv[1]) != 3:
    d = "" 
    argv = []  
    while len(d) != 1 and len(d) != 3:
        d = str(raw_input("Something is wrong. The program can not use your data. Please try to write them in other way."))
    argv = argv + [""] +[d]
                                       
if len(argv) != 1:
    amc = argv[1]
    if len(amc) == 1:
        a = amc.upper() 
        a = func(dicaaa, a)
        if a != 0:
            print a
        if a == 0:
            print "Xxx"     
    if len(amc) == 3:
        b = amc[0].upper() +  amc[1:3]
        b = func(dica, b)
        if b != 0:
            print b
        if b == 0:
            print "X"  
    