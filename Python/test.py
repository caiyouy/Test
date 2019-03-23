#! /usr/bin/python
import math
def func(x):
    y=x**3.0+2.0/3.0*x**2.0+1.0/9.0*x-4.0/9.0
    return -y

x=(math.sqrt(17)-1.0)/6.0
y=func(x)
print x
print y
