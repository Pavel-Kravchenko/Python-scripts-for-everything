import random
print "Happy new year!!!!!"
print "Write your name please"
n = raw_input() 
print "Here you can tipe your congratulations"
a = str(raw_input())
n = str(n)[0].upper() + str(n)[1:]
r = open ("congratulations.txt", "r")
ty = r.readline()
w = ""
while len (ty) >0:
    w = w+ ty
    ty = r.readline()
 

r.close() 
q = open("congratulations.txt", "w")
q.write(w)
q.write(n + "\t" + "wrote:" + "\n")
q.write("'"+a +"'"+ "\n")
print "Thank you for your words"
print "Have a great time!"
price  = ["!","@","$","^","&","*",":",")","O",">",]
r = len(a) * 100
e=""
while r !=0 :
    f = random.choice(price) 
    e = e + f
    r -= 1
print e
q.close()