
# coding: utf-8

# In[66]:

input_n = int(input().strip())
YES = {"-"}
NO = {"-"}
input_string_beatrice_1 = ""
input_string_beatrice_2 = ""
if input_n < 1000:
    while input_string_beatrice_1 != "HELP":
        input_string_beatrice_1 = input().strip()
        if input_string_beatrice_1 == "HELP":
            break
        input_string_beatrice_2 = input().strip()
        if input_string_beatrice_2 == "NO":
            _set_ = set(input_string_beatrice_1.split(" "))
            NO = NO.union(_set_)
        if input_string_beatrice_2 == "YES":
            _set_ = set(input_string_beatrice_1.split(" "))
            YES = YES.union(_set_)
            YES = _set_.intersection(YES)
_out_ = YES.difference(NO)
_out_ = sorted(list(_out_))
print(" ".join(_out_))

