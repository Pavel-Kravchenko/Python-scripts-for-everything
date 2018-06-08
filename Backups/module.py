def worker(k):

    i=0                                                                   
    j=3                                                                   
    triplets = []                                                         
    while j < len(k)+1:    
       # here we have a controvercial var                                                 
        triplets = triplets + [k[i:j]]                                    
        i = i + 3                                                         
        j = j + 3                                                                                                                              
    G = 0                                                                 
    C = 0                                                                 
    A = 0                                                                 
    T = 0                                                                 
    for w in triplets:                                                    
        if w[2] == "G":                                                   
            G = G + 1                                                     
        if w[2] == "C":                                                   
            C = C + 1                                                     
        if w[2] == "T":                                                   
            T = T + 1                                                     
        if w[2] == "A":                                                   
            A = A + 1                                                     
                                                                          
    sa = A + T + G + C                                                   
    al = G + C  
    if al == 0:
        return(0)                                                         
    r = float(al)/float(sa) * 100                                       
    r = round(r)                                                          
    r = int(r)                                                            
    return(r)                                                         
    
