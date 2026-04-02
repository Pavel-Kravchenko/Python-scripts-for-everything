import random
def i(number):
    nuc = ["A", "T", "G", "C"]
    calc = 0
    line = ""
    while calc < number:
        lis = str(random.choice(nuc)) 
        line = line + lis
        calc = calc + 1
    return line

def ii(num):
    import Kravchenko_pr9_randomdnamany
    line = Kravchenko_pr9_randomdnamany.tt(num)
    X = 0                       
    Y = 60
    q = ""                      
    if len(line) > 60:          
        d = len(line)           
        while d > 0:            
            j = line[X:Y]       
            d = d - 60          
            if d > 0:
                q = q + str(j) + "\n" 
   
            else:
                q = q + j
           
            X = X + 60          
            Y = Y + 60          
    else:
        q = line                       
    if num == 0:
        q = ">random-" + str(num + 1) + "\n" + str(q)
    else:
        q = "\n" + ">random-" + str(num + 1) + "\n" + str(q)   
    return q