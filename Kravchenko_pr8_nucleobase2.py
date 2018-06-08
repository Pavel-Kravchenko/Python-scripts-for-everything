inp = raw_input()

pyrimidine = ["T", "C" ,"U"]
purine = ["A", "G"]
inp = inp.upper()
if inp in pyrimidine :
    print ("pyrimidine")
if inp in purine :
    print ("purine")
else :
    print ("invalid data")


