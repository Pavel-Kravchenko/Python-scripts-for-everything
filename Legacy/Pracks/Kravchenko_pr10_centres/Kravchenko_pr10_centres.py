from sys import argv
indata = open(argv[1], "r")
line = indata.readline()
line = indata.readline()
line = line.strip()
import erth
kol= 1
pu = ""
fila = open("basic.tab", "w")
pol = []
summ = 0
summm = 0
summmm = 0
while len(line) != 0:
    l = line[4]
    if l not in ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "V", "U", "W", "X", "Y", "Z"]:
        print "invalid data" 
        raise SystemExit(1)
    if kol > 0:
        if l != pu:
            pu = l                                                                                                              
            eda = summ / kol                                              
            eda2 = summm / kol                                            
            eda3 = summmm / kol                                           
            rty = str(eda) + "\t" + str(eda2) + "\t" + str(eda3) + "\n"   
            fila.write(rty)                                               
            summ = 0                                                      
            summm = 0                                                     
            summmm = 0                                                    
            kol = 0                                                                                                           
    pol = erth.calc(line)
    a = pol[0]
    b = pol[1]
    c = pol[2] 
    summ = summ + a                
    summm = summm + b
    summmm = summmm + c
    line = indata.readline()         
    line  = line.strip()              
    kol = kol + 1
eda = summ / kol                                              
eda2 = summm / kol                                            
eda3 = summmm / kol                                           
rty = str(eda) + "\t" + str(eda2) + "\t" + str(eda3) + "\n"   
fila.write(rty)   
fila.close()                    
indata = open("basic.tab", "r")
line = indata.readline()
line = indata.readline()
while len(line) > 0: 
    print line
    line = indata.readline()
   


