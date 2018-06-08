def func(l):
    y = l         
    b = ""              
    for i in l:
        if i == "#":
            break     
        if i != "#":
            b = b + i
    return b
  