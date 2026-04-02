inp = raw_input("Write your amino acid. Example: CH-H     ")
inp = inp.upper()
print inp
tiny = ["G", "A", "C", "T", "S", "CH-H", "S"]
aliphatic = ["V", "I", "L"]
small = ["CS-S", "P", "V", "D", "N", "G", "S", "A", "CH-H", "T", "D"]
positive = ["K", "R", "H"]
polar = ["Q", "S", "G", "CH-H","T", "D", "Y", "W", "H", "K" ,"R", "E"]
charged = ["D","K", "R" ,"H", "E"]
non_polar = ["N", "CS-S", "G", "A", "CH-H", "V", "I", "L", "M", "F", "Y", "W", "H", "K"]
aromatic = ["F", "Y", "W", "H"]
if inp in tiny or aliphatic or small or positive or polar or charged or non_polar or aromatic:
    if inp in tiny:
        print "tiny"
    if inp in aliphatic:
        print "aliphatic"
    if inp in small:
        print "small"
    if inp in positive:
        print "positive"
    if inp in polar:
        print "polar"
    if inp in charged:
        print "charged"
    if inp in non_polar:
        print "non-polar"
    if inp in aromatic:
        print "aromatic" 



