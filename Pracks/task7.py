
# coding: utf-8

# In[ ]:

time = int(input())
line = input()
lenl = len(line)
def liner(d):
    dic = {}
    for k in d:
        if k >= 0 and k < lenl: 
            if k == 0:
                if d[k] == "B" and d[k+1] == "G":
                    dic[k] = "G"
                if d[k] == "B" and d[k+1] == "B":
                    dic[k] = "B"
                if d[k] == "G":
                    dic[k] = "G"
            if k > 0 and k < lenl-1:
                if d[k] == "B" and d[k+1] == "G":
                    dic[k] = "G"
                if d[k] == "B" and d[k+1] == "B":
                    dic[k] = "B"
                if d[k] == "G" and d[k-1] == "B":
                    dic[k] = "B"
                if d[k] == "G" and d[k-1] == "G":
                    dic[k] = "G"
            
            if k+1 == lenl:
                if d[k] == "G" and d[k-1] == "B":
                    dic[k] = "B"
                if d[k] == "G" and d[k-1] == "G":
                    dic[k] = "G"
                if d[k] == "B":
                    dic[k] = "B"         
    d = {}
    for key in sorted(dic.keys()):
        d[key] = dic[key]
    return(d)

if time <= 50 and lenl > 1:
    d = {}
    for i in range(len(line)):
        d[i] = line[i]
    c = 0
    while c != time:
        d = liner(d)  
        c += 1
        
    line2 = ""
    for k in d:
        line2 = line2 + d[k]
    print(line2)
if time <= 50 and lenl == 1:
    if line == "B":
        print("B")  
    if line == "G":
        print("G")  

