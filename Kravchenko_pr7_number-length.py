indata = open("numbers.txt", "r")
s1 = indata.readline()
s2 = indata.readline()
b =  int(s1) ** int(s2)
a = str(b)
print len(a)
indata.close()