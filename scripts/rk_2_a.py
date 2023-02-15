import numpy as np
import math as m
from tabulate import tabulate

"""
DUE NEXT WEDNESDAY 02/22/2023
RK_2_A 
RK_2_B 
RK_4_A 
RK_4_B 
"""

"""
R_K_2_B AND C 
Format short g 
T1 R1  R2  R3 
0.2  0.2618 0.25574 0.25282 
2.0  0.26105 0.25495 0.25232 
"""

"""
R_K_4_B output 
T1   R1     R2    R3
0.2  15.33  15.67 15.836
2.0  15.586 15.81 15.91
"""

actual_function = lambda x: (x+1)*(x+1) - 0.5*m.exp(x)
function = lambda x, y: y - x**2 + 1

def RK2(a:float, b:float, N:int, x0: float, y0: float):
    """Computes RK2 for a given function"""
    x = np.zeros(N+1)
    y = np.zeros(N+1)
    actual_vals = np.zeros(N+1)
    x[0] = x0
    y[0] = y0
    actual_vals[0] = actual_function(x0)
    h = (b-a)/N

    for i in range(N):
        x[i+1] = x[i] + h
        s1 = function(x[i], y[i])
        s2 = function(x[i+1], y[i]+(h*s1))
        y[i+1] = y[i] + (h*(s1 + s2)/2)
        actual_vals[i+1] = actual_function(x[i+1])
    return x, y, actual_vals

def RK4(a:float, b:float, N:int, x0:float, y0:float, actual_function=function,
    function=function):
    """Computes RK4 for a given function"""
    x = np.zeros(N+1)
    y = np.zeros(N+1)
    actual_vals = np.zeros(N+1)
    x[0] = x0
    y[0] = y0
    actual_vals[0] = actual_function(x0)
    
    #round to ceiling of h
    h = (b-a)/N
    print("h = ", h)

    for i in range(N):
        x[i+1] = x[i] + h
        s1 = function(x[i], y[i])
        s2 = function(x[i] + h/2, y[i] + (h*s1)/2)
        s3 = function(x[i] + h/2, y[i] + (h*s2)/2)
        s4 = function(x[i] + h, y[i] + h*s3)
        y[i+1] = y[i] + (h*(s1 + 2*s2 + 2*s3 + s4)/6)
        actual_vals[i+1] = actual_function(x[i+1])

    return x, y, actual_vals

def R_K2_A():
    """RK2A Due Next Wednesday """
    a = 0 
    b = 2
    ad = 0.5
    N = 10

    x, y, actual_vals = RK2(a, b, N, a, ad)
    print("RK2A")
    for i in range(N+1):
        print("x = ", x[i], "y = ", y[i], "actual = ", actual_vals[i])

def RK4_A():
    """RK4A Due Next Wednesday """
    a = 0 
    b = 2
    ad = 0.5
    N = 10

    x, y, actual_vals = RK4(a, b, N, a, ad)
    print("RK4A")

    #print in table format
    for i in range(N+1):
        print("x = ", x[i], "y = ", y[i], "actual = ", actual_vals[i])


R_K2_A()
RK4_A()


## RK45C
a = 0
b = 2
ad = 0
N = 10
k = 6.22E-19

    

