indata = open("cds.txt", "r")
s1 = indata.readline()
s2 = indata.readline()
m = min(int(s1),int(s2))
n = max(int(s1),int(s2))
print (n - m - 2)/3 
indata.close()