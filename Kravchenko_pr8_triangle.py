import math 
a = float(raw_input())
b = float(raw_input())
c = float(raw_input())

x = max(a, b, c)

z = (-b ** 2 - c ** 2 + a ** 2)/(-1*(2 * b * c))
e = (-a ** 2 - c ** 2 + b ** 2)/(-1*(2 * a * c))
c = (-b ** 2 - a ** 2 + c ** 2)/(-1*(2 * a * b))

u =  math.acos(z)
h =  math.acos(e)
n =  math.acos(c)
"""if x == (a**2 + b**2)**0.5:   
    print "right" """            
if x == a + b:
    print "degenerate"
    
if x > a + b:
    print "impossible"
o = max(u, h, n)
d = math.pi
if x < a + b:
    if o > d/2:
        print "obtuse"
    if o < d/2:
        print "acute"
    if o == d/2:
       print "right"




