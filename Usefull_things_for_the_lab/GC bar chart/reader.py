from sys import argv


# collecting data in reding module and the table making
name_of_1 = argv[1]
in_file_1 = open(name_of_1, "r")

name_of_out_file = argv[2]
out_file = open(name_of_out_file, "w")

out_file.write("Start" + "\t" + "End" + "\t"+ "Name" + "\n") 

# start readin and writing files                             
# module
#///////////////////////////////////////////////////////////////////
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
       
        if w[0:2]=="GC" or w[0:2]=="CT" or w[0:2]=="AT" or w[0:2]=="GT" or w[0:2]=="TC" or w[0:2]=="CC" or w[0:2]=="AC" or w[0:2]=="GG" or w[0:2]=="CG":                                                    
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
    print sa, al                                                        
    r = float(al)/float(sa) * 100                                       
    r = round(r)                                                          
    r = int(r)
    print r                                                            
    return(r)                                                         
    
#/////////////////////////////////////


def GC_calculator(p):
    radius = 3
    line = []                                                                       
    while len(p) > 0:                                                                      
        line = line + [p[:radius]]                                                  
        p = p[radius:]                                                              
                                                                                    
    # some check of elements                                                        
    new_line = []                                                                   
    for k in line:                                                                  
        if len(k) == radius:                                                        
            new_line = new_line + [k]                                               
        if len(k) < radius:                                                         
            while len(k) > (radius - 3):                                            
                k = k[:-1]                                                          
            new_line = new_line + [k]                                               
    # restriction of 0- elements                                                    
    new_line1 = []                                                                  
    for k in new_line:                                                              
        if len(k) != 0:                                                             
            new_line1 = new_line1 + [k]                                             
                                                                                    
    # we have ready sequence now!                                                   
    # let's use a module for calculations                                           
                                                                                    
    dic = {}                                                                        
    calc = 0 
    ou = 0 
                                                                      
    for k in new_line1:
        out_module = worker(k)                                                      
        if out_module != 0:
            ou = ou + out_module
            calc = calc + 1
    print ou, calc
    if ou != 0:
        dic = float(ou)/float(calc)
    else:
        dic = 0
    print dic    
    return dic
                                                                                
#///////////////////////////////////////////////////////////
                                                                                
#///////////////////////////////////////////////////////////
                                                                                 
s = in_file_1.readline()                                                  
s = s.strip()
s = s.split(" ")                                                                                                                             
name_of_file = s[7]
s = in_file_1.readline()                                                  
s = s.strip()
ww = {}
colc = 1
# start the bacic sycle and making the table of seqences 
while len(s) != 0:
    s = s.split(" ")   
    # lookind for "gene" in line                                     
    if s[0] == "gene" and s[1] == "":
        length = s[12] 
        # in cause of complement using this condition
        if "complement" in length:
            length = length.split("complement")
            length = length[1].split("..") 
            w = ((length[0])[1:])
            e = ((length[1])[:-1])
                    
        if "complement" not in length and "(" not in length[0] :
            length = length.split("..")                        
            w = (length[0]) 
            e = (length[1])
        
        s = in_file_1.readline()                                              
        s = s.strip()
        # looking for sequence name
        if "/locus_tag=" in s:                                    
            name_of_gene = s[12:-1]
        
        if "/gene=" in s:  
            s = in_file_1.readline()                                              
            s = s.strip()
            name_of_gene = s[12:-1]
        e = str(e) + "\t" + name_of_gene + "\t" + str(colc)
        colc = colc + 1
        
        if w not in ww.keys():
            ww[w] = e 
    s = in_file_1.readline()                                              
    s = s.strip() 
ww = sorted(ww.keys())
print ww

for i in ww:
    out_file.write(i + "\t" + ww[i] + "\n")                                                        
out_file.close()
in_file_1.close()


# module for collecting sequence
#////////////////////////////////////////////////////////////
name_of_1 = argv[1]
in_file_1 = open(name_of_1, "r")
s = in_file_1.readline()                                                  
s = s.strip()
while len(s) != 0:
    if "ORIGIN" in s:                 
        p = "" 
        s = in_file_1.readline()
        s = s.strip()
        while "//" not in s:
            s = s.split(" ")
            a = s.pop(0)
            if len(s) == 1:
                p = p + s[0]
            if len(s) == 2:
                p = p + s[0] + s[1]
            if len(s) == 3:
                p = p + s[0] + s[1] + s[2]
            if len(s) == 4:
                p = p + s[0] + s[1] + s[2] + s[3]
            if len(s) == 5:
                p = p + s[0] + s[1] + s[2] + s[3] + s[4]
            if len(s) == 6:
                p = p + s[0] + s[1] + s[2] + s[3] + s[4] + s[5]
            s = in_file_1.readline()
            s = s.strip()
    else:
        s = in_file_1.readline()
        s = s.strip()
p = p.upper()


# starting a new module
#////////////////////////////////////////////////////////////////////////
# it will find a sequence of genes and write them in other dictionary
in_file_2 = open(name_of_out_file, "r")
di = {}
de = {}
y = 1

s = in_file_2.readline()                                                  
s = in_file_2.readline()                                                  
s = s.strip()
split = s.split("	")  
r = int(split[0]) - 1              
f = int(split[1])  

if r != 0:                  
    left = 0               
    right = r  
    srez = p[left:right]
    de[y] = srez
    y = y + 1
    # making two dict for introns and exones  "di" and "de"     
    while len(s) != 0:
        l = p[r:f]              
        di[y] = l 
        y = y + 1               
        s = in_file_2.readline() 
        if len(s) == 0:
            break                                                 
        else:
            s = s.strip()
            split = s.split("	")     
        left = f                
        right = int(split[0]) - 2  
        srez = p[left:right]                                 
        de[y] = srez 
        y = y + 1
        r = int(split[0]) - 1              
        f = int(split[1])
    l = p[r:f]              
    di[y] = l       

if r == 0:                  
    l = p[r:f]
    di[y] = l 
    y = y + 1
    # making two dict for introns and exones  "di" and "de"     
    while len(s) != 0:
        s = in_file_2.readline() 
        if len(s) == 0:
            break                                                 
        else:
            s = s.strip()
            split = s.split("	")     
        left = f                
        right = int(split[0]) - 1   
        srez = p[left:right]
        de[y] = srez 
        y = y + 1
        r = int(split[0]) - 1              
        f = int(split[1])
        l = p[r:f]
        di[y] = l 
        y = y + 1
    l = p[r:f]              
    di[y] = l       
 
#//////////////////////////////////////////////////    
out_file1 = open("Out_sequence_in_intrones_and_exones.txt","w")
out_file2 = open("Out_DC_percent.txt","w")
h = int(len(de) + len(di))
calc = 0
calce = 0
calci = 0
for i in range(h):  
    if i in di:                    
        dic = GC_calculator(di[i])
        calci = calci + 1  

        out_file1.write(str(i))       
        out_file1.write("\t" + di[i] + "\n") 

in_file_2.close()                                 
in_file_1.close()
print "Check your file"        
