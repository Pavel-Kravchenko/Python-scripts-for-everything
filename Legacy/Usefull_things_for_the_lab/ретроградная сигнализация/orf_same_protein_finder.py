from sys import argv

name_of_1 = argv[1]
name_of_2 = argv[2]
name_of_out_file = argv[3]

in_file_1 = open(name_of_1, "r")
in_file_2 = open(name_of_2, "r")
out_file = open(name_of_out_file, "w")

line1 = in_file_1.readline()
line1 = line1.strip()

line2 = in_file_2.readline()
line2 = line2.strip()

a=[]
b=[]
                               
while len(line1) > 0:
    if line1[0] != ">":
       a = a + [line1]
    line1 = in_file_1.readline()
    line1 = line1.strip()

while len(line2) > 0:
    if line2[0] != ">":
       b = b + [line2]
    line2 = in_file_2.readline()
    line2 = line2.strip()

calc = 0

for i in a:
    if i in b:
        calc = calc + 1
        out_file.write(i + "\n") 

print "Check your file"        
print "There have been",calc, "objects found"
in_file_1.close()
out_file.close() 