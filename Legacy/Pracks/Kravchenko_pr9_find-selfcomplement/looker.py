def look(name, op):                                           
    indata = open(name, "r")                                               
    line = indata.readline()                   
    qwerty = ""                      
    while len(line) > 0:                       
        line = line.strip()                    
        if len(line) > 0:                      
            if str(op) in line:
                             
                while len(line) > 0:
                    line = indata.readline()
                    line = line.strip()
                    if ">" in line:
                        break
                    qwerty = qwerty + line
                       
        line = indata.readline()               
    indata.close()
    return qwerty                             
                                               