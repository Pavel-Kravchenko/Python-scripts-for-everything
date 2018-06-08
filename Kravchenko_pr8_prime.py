a = int(raw_input())
b = a%2 
c = a%3
d = a%4
e = a%5
f = a%6
g = a%7
h = a%8
i = a%9 
if a == 3 or a == 5 or a == 7 or a == 2:
    print "Prime"             
else:              
    if a == 0 or a == 1 or b == 0 or c == 0 or d == 0 or e == 0 or f == 0 or g == 0 or h == 0 or i == 0:
        print "Composite"
    else:
        print "Prime"  
