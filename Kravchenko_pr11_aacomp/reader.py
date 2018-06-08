def func(data):
    d = {}
    a = open(data,"r")
    line = a.readline()
    d = {}
    for e in a:
        for i in e:    
            if i.isalpha() == True:          
                i = i.upper()       
                if i not in d:      
                    d[i] = 1        
                else:               
                    d[i] += 1                   
    return d
