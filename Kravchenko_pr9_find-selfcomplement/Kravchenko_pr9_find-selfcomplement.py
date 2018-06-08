from sys import argv
a = str(argv[1])
b = str(argv[2])
c = str(argv[3])
fila = open(c, "w")

import looker
assa = looker.look(a,b)
import palindromcalc

fila.write(">")
fila.write(b + "\n")

da = palindromcalc.pal(assa)
for i in da:
    fila.write(i + "\n")