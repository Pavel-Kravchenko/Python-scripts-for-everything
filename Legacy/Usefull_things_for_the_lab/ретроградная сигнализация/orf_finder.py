from sys import argv

if argv[1] == None:
    name_of_in_file = raw_input("Enter the name of in file")             
else: 
    name_of_in_file = argv[1]
    
if argv[2] == None:
    name_of_out_file = raw_input("Enter the name of out file")
else: 
    name_of_out_file = argv[2]

min_a = raw_input("Enter the min aminoacids number")
max_a = raw_input("Enter the max aminoacids number")

in_file = open(name_of_in_file, "r")
out_file = open(name_of_out_file, "w")

line = in_file.readline()
line = line.strip()
a=[]
                               
while len(line) > 0:
    if line[0] != ">":
       a = a + [line]
    line = in_file.readline()
    line = line.strip()

c = 0

for i in a:
    if len(i) > int(min_a) and len(i) < int(max_a):
        out_file.write(i + "\n") 
        c = c + 1

print "Check your file. There have been ",c," objects found"        
in_file.close()
out_file.close() 