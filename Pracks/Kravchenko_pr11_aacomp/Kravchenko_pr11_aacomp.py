from sys import argv
import reader
if len(argv) == 1:
    a = ""
    while len(a) == 0:
        a = str(raw_input("Write the name of in-file please"))
if len(argv) > 1:
    if len(argv[1]) > 0:       
        a = str(argv[1])
dic = reader.func(a)
d = sorted(dic.keys()) 

if len(argv) > 1:  
    if len(argv) == 3:
        b = str(argv[2]) 
        b = open(b,"w")
        for i in d:
            b.write(str(i) + "\t" + str(dic[i]) + "\n")
        b.close()                
    else:
        print "The total is:"
        for i in d:
            print str(i) + "\t" + str(dic[i]) 
    raise SystemExit(1)  
               
print "The total is:"                  
for i in d:                            
    print str(i) + "\t" + str(dic[i])  