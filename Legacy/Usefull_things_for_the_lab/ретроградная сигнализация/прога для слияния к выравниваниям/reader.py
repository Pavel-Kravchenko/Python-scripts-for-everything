from sys import argv

name_of_1 = argv[1]
name_of_out_file = argv[2]

in_file_1 = open(name_of_1, "r")
out_file = open(name_of_out_file +".fasta", "w")

s = in_file_1.readline()
s = in_file_1.readline()
p = ""
while len(s) != 0:
    p = p + s
    s = in_file_1.readline()
    

import win32com.client
Excel = win32com.client.Dispatch("Excel.Application")

wb = Excel.Workbooks.Open(u'C:\\Users\\Andrey\\Desktop\\Павел\\Биоинформатика\\ретроградная сигнализация\\прога для слиянию к выравниваниям\\xl.xls')
sheet = wb.ActiveSheet
val = 1
sheet.Cells(1,2).value = val
i = 1
for rec in vals:
    sheet.Cells(i,3).value = rec
    i = i + 1
wb.Save()
wb.Close()
Excel.Quit()

out_file.write(s)
in_file_1.close()
print "Check your file"        
out_file.close() 
