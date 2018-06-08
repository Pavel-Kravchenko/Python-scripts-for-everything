def pal(assa):
    a = []
    for i in range(len(assa)):
        if i > 0 & i < len(assa): 
            left = len(assa[:i]) - 1               
            right = len(assa[:i]) + 1
            if len(assa) > right:                                    
                if assa[left] == assa[right] or len(assa) > right:               
                    while assa[left] == assa[right]:        
                        if len(assa[left:right + 1]) > 3:   
                            a = a + [assa[left:right + 1]]  
                                                            
                            print a                         
                        left = left - 1                     
                        right = right + 1
                        if len(assa) == right:
                            break                   
                        print a                             
                                                            
    return a
    
