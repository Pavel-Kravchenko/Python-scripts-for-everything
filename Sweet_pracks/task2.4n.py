
# coding: utf-8

# In[14]:

input_n = int(input().strip())
_set_ = {"--"}
for i in range(1, input_n + 1):
    _set_ = _set_.union(set(str(i)))
input_string_beatrice_1 = ""
input_string_beatrice_2 = ""
if input_n < 1000:
    while input_string_beatrice_1 != "HELP":
        input_string_beatrice_1 = input().strip()
        if input_string_beatrice_1 == "HELP":
            break
        input_string_beatrice_2 = input().strip()
        if input_string_beatrice_2 == "NO":
            set_NO = set(input_string_beatrice_1.split(" "))
            _set_ = _set_.difference(set_NO)
        if input_string_beatrice_2 == "YES":
            set_YES = set(input_string_beatrice_1.split(" "))
            _set_ = _set_ & set_YES
_out_ = sorted(list(_set_))
_out_e = " ".join(_out_)
print(str(_out_e))

