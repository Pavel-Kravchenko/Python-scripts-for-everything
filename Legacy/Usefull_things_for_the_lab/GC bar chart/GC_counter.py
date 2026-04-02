from sys import argv

# open everything
name_of_1 = argv[1]
name_of_out_file = argv[2]
radius = int(raw_input("Type your radius in amino acids, please: "))
radius = radius*3

in_file_1 = open(name_of_1, "r")
out_file = open(name_of_out_file +".txt", "w")

# retrieving data                                                                   
s = in_file_1.readline()                                                  
s = in_file_1.readline()                                                  
s = s.strip()                                                             
p = ""                                                                    

while len(s) != 0:                                                        
    p = p + s                                                             
    s = in_file_1.readline()                                              
    s = s.strip()                                                         
# maiking a upper unique! line
p = p.upper()
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
import module

dic = {}
calc = ""

for k in new_line1:
    out_module = module.worker(k)
    dic[len(calc)] = out_module
    calc = calc + str(1)

for i in dic:
    out_file.write(str(i) + "\t" + str(dic[i]) + "\n")
in_file_1.close()
print "Check your file"        
 
