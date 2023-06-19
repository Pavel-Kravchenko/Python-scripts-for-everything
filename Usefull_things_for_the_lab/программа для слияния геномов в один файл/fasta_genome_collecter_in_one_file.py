from sys import argv

name_of_1 = argv[1]
name_of_out_file = raw_input("Type the name of out file: ")

def module_writer(mon):
    X = 0                                 
    Y = 60                            
    if len(mon) > 60:                   
        d = len(mon)                    
        while d > 0:                  
            j = mon[X:Y]                
            d = d - 60                
            out_file.write(j + "\n" ) 
            X = X + 60                
            Y = Y + 60                
    else:                             
        out_file.write(mon + "\n" )     

   
in_file_1 = open(name_of_1, "r")
s = in_file_1.read()

h = 0
for i in s:
    if i == ">":
        h = h + 1 
    
in_file_1.close()
in_file_1 = open(name_of_1, "r")

line1 = in_file_1.readline()        
line1 = line1.strip()

import os
name_of_out_file = name_of_out_directory
os.mkdir(name_of_out_directory)
os.chdir(name_of_out_file)
         
calc = 1               
while len(linel) != 0:
    out_file = open(name_of_out_file + str(calc) +".fasta", "w")
    a=""
    if line1[0] == ">": 
        out_file.write(line1 + "\n") 
        line1 = in_file_1.readline()        
        line1 = line1.strip() 
 
    if line1[0] != ">" and line1[0] != 0: 
        while line1[0] != ">":
            a = a + line1                     
            line1 = in_file_1.readline()        
            line1 = line1.strip() 
            if len(line1) == 0:
                break    
    module_writer(a)                 

    b="" 
    if line2[0] == ">": 
        out_file.write(line2 + "\n") 
        line2 = in_file_2.readline()   
        line2 = line2.strip() 
         
    if line2[0] != ">" and line2[0] != 0: 
        while line2[0] != ">":                       
            b = b + line2                    
            line2 = in_file_2.readline()         
            line2 = line2.strip()
            if len(line2) == 0:
                break            
    module_writer(b)
    out_file.close()
    calc = calc + 1 
    h = h - 1              
                                                
print "Check your file"        
in_file_1.close()
in_file_2.close()
out_file.close() 