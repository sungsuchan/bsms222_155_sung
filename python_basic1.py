# -*- coding: utf-8 -*-
"""
practice python in spyder.

"""

age = 22
name = 'suchan'
#%%
print('{0} was {1} years old when he worte this book'.format(name,age))
print('Why is {0} playing with that python?'.format(name))
#%%
print(f'{name} was {age} years old when he wrote this book')
print(f'why is {name} playing with that python?')

#%%
print('{0:.3f}'.format(1.0/3))
print('{0:_^11}'.format('hello'))
#%%
print('a', end='')
print('b', end='')
#%% 
print('a', end=' ')
print('b', end =' ')
print('c', end=' ')
#%%
print('this is the first line \nthis is the second line')
print('this is the first line \t this is the second line')
#%%
print(R'new \n not line lol')
#%%
i=5
print(i)
i = i+1
print(i)

s= '''This is a multi-line string.
This is the second line.'''
print(s)
#%%
i=5
print(i)

i=5; print(i)
#%%
s='this is a string. \
this continues the string.'
print(s)
#%%
i = \
5
print(i)
#%%
'la' * 3
13 % 3 #returns the remainder of the devision.
#%%
a =2 
a *= 3
#%%
length = 5
breadth = 2

area = length * breadth
print('Area is', area)
print('Perimeter is', 2*(length + breadth))









































