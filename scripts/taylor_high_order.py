import numpy as np


h = 0.2 
y0 = 0.5 
x0 = 0.0
N = 10 

def function(x, y):
    return y - x**2 + 1

def derivative(x, y):
    return y - x**2 - 2*x + 1


for i in range(N):
    k1 = function(x0, y0)
    k2 = derivative(x0, y0)

    y0 = y0 + h*k1 + h**2/2*k2
    
    x0 = x0 + h
    print("x = ", x0, "y = ", y0)