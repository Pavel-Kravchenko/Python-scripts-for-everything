from sys import argv

name_of_1 = argv[1]
name_of_out_file = argv[2]

in_file_1 = open(name_of_1, "r")
out_file = open(name_of_out_file +".txt", "w")

s = in_file_1.readline()
s = in_file_1.readline()
s = s.strip()
p = ""
while len(s) != 0:
    p = p + s
    s = in_file_1.readline()
    s = s.strip() 
out_file.write(p)
in_file_1.close()
print "Check your file"        
out_file.close() 
