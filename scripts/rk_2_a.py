import numpy as np
import math as m

function = lambda x, y: y - x**2 + 1

def RK2(a:float, b:float, N:int, x0: float, y0: float):
    """Computes RK2 for a given function"""
    x = np.zeros(N+1)
    y = np.zeros(N+1)
    x[0] = x0
    y[0] = y0
    h = (b-a)/N
    for i in range(N):
        x[i+1] = x[i] + h
        s1 = function(x[i], y[i])
        s2 = function(x[i+1], y[i]+(h*s1))
        y[i+1] = y[i] + (h*(s1 + s2)/2)
    return x, y

def RK4(a:float, b:float, N:int, x0:float, y0:float):
    """Computes RK4 for a given function"""
    x = np.zeros(N+1)
    y = np.zeros(N+1)
    x[0] = x0
    y[0] = y0
    
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

    return x, y

def R_K2_A():
    """RK2A Due Next Wednesday """
    a = 0 
    b = 2
    ad = 0.5
    N = 10

    x, y = RK2(a, b, N, a, ad)
    print("RK2A")
    for i in range(N+1):
        print("x = ", x[i], "y = ", y[i])

def RK4_A():
    """RK4A Due Next Wednesday """
    a = 0 
    b = 2
    ad = 0.5
    N = 10

    x, y = RK4(a, b, N, a, ad)
    print("RK4A")
    for i in range(N+1):
        print("x = ", x[i], "y = ", y[i])
        
R_K2_A()
RK4_A()

"""
R_K_2_B AND C - Invoke Euler 
Format short g 
T1 R1  R2  R3 
0.2  0.2618 0.25574 0.25282 
2.0  0.26105 0.25495 0.25232 
"""

    

