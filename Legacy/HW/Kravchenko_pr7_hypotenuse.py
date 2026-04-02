indata = open("catheti.txt", "r")
s1 = indata.readline()
s2 = indata.readline()
a = float(s1)
b = float(s2)
h = a ** 2 + b ** 2
f = 0.5
c = h ** f
print c 
indata.close()