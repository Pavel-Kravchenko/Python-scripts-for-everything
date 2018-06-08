def func(data):
    d = {}
    a = open(data,"r")
    line = a.readline()
    d = {}
    for e in a:
        if e[0] != ">":
            for i in e:                    
                if i.isalpha() == True:      
                    i = i.upper()          
                    if i not in d:         
                        d[i] = 1           
                    else:                  
                        d[i] += 1    
        else:
            name = srt(e[1:]) + ".tab"
            pond = open(name, "w")
            pond.write(e + "\n")
            pond.write(e + "\n")
            p = sorted(d.keys())   
    return d
