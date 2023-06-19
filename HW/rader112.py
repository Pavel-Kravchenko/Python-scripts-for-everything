from sys import argv


# collecting data in reding module and the table making
name_of_1 = argv[1]
in_file_1 = open(name_of_1, "r")

name_of_out_file = argv[2]


out_file2 = open(str(2) + name_of_out_file, "w")
 
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
         
       # d = {}
       # y = 1
       # for i in p:
        #    d = dict([(y, i)])
         #   y = y + 1
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
 
        
out_file = open("txt.txt","w")
h = int(len(de) + len(di))
for i in range(h):
    if i in de:
        out_file.write("introne" + "\t")       
        out_file.write(str(i))       
        out_file.write("\t" + de[i] + "\n")       
    if i in di:                    
        out_file.write("exone" + "\t")       
        out_file.write(str(i))       
        out_file.write("\t" + di[i] + "\n") 

      

in_file_2.close()                                 
in_file_1.close()
print "Check your file"        
