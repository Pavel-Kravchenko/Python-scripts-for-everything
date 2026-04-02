def wolf(data):
    data = str(data)
    data = data.strip()
    if data[0:6] == "ATOM  " and data[12:16] == " CA ":
        dat = data[17:55]       
        dat = dat. split()          
        dat = "\t ".join(dat)                                         
    if data[0:6] != "ATOM  " or data[12:16] != " CA ":
        dat = ""
    return dat

                            