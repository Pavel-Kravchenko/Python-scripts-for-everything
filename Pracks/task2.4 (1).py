
# coding: utf-8

# In[47]:

input_n = int(input().strip())
YES = set()
YES2 = set()
NO = set()
NO2 = set()
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
            if len(NO) != 0:
                NO = NO.intersection(_set_)
                print(NO)
            if len(NO) == 0:
                NO = NO.union(_set_)
                print(NO)
        if input_string_beatrice_2 == "YES":
            _set_ = set(input_string_beatrice_1.split(" "))
            if len(YES) != 0:
                YES = YES.intersection(_set_)
                print(YES)
            if len(YES) == 0:
                YES = YES.union(_set_)
                print(YES)
_out_ = YES.difference(NO)
_out_ = sorted(list(_out_))
print(" ".join(_out_))

