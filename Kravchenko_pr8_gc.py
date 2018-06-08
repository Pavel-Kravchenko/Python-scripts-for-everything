from sys import argv
indata = open(argv[1], "r")
line = indata.readline()
line = line.strip()
a = ""
while ">" in line:
    line = indata.readline()
    line = line.strip()                                  
while len(line) > 0:
    a = a + line
    line = indata.readline()
    line = line.strip()
a = a.upper()
s = 0
d = []
y = 0
h = 0
k = 0
m = 0

p = a
while len(p) != s:
    e = a[s]
    e = [e]
    d = d + e  
    s = s + 1
while "A" in d:
    for i in d:
        if i == "A":
            y = y + 1
            d.remove(i)
while "T" in d:
    for i in d:
        if i == "T":
            h = h + 1
            d.remove(i)
while "G" in d:
    for i in d:
        if i == "G":
            k = k + 1
            d.remove(i)
while "C" in d:
    for i in d:
        if i == "C":
            m = m + 1
            d.remove(i)
sas = y + h + k + m
aal = k + m
r = float(aal)/float(sas) * 100

r = round(r)

r = int(r)
print r, "%"

indata.close()