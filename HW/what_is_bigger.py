indata1 = open("numbers1.txt", "r")

s1 = indata1.readline()
s2 = indata1.readline()
indata1.close()
indata2 = open("numbers2.txt", "r")

s11 = indata2.readline()
s21 = indata2.readline()
indata2.close()

a =  int(s1) ** int(s2)
b =  int(s11) ** int(s21)
t = " The first one ="
m = "The second one ="
print t,  a, m, b  
z = max(a,b)

w= "The bigger one is"
print w, z
