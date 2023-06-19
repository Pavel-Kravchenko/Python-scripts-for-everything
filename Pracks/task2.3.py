
# coding: utf-8

# In[ ]:

_input_ = str(input()).strip()
_input_ = _input_.split(" ")
_list_ = ''

for i in range(int(_input_[0])):
    _list_ = _list_ + "I"

for i in range(int(_input_[1])):
    block = input().strip().split()
    left = int(block[0])
    right = int(block[1])
    if left >= 1 and right <= 100 and left <= 100 and right >= 1:
        if left <= len(_list_) and right <= len(_list_):
            if left - right != 0:
                for i in range(left, right + 1):
                    _list_ = _list_[:(i - 1)] + "." + _list_[(i):]
            else:
                _list_ = _list_[:(left - 1)] + "." + _list_[(right):]
print(_list_)

