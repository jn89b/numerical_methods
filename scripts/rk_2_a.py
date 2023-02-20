import numpy as np
import math as m
from tabulate import tabulate

import taylor_high_order


"""
DUE NEXT WEDNESDAY 02/22/2023
RK_2_A - HW 7
RK_2_B - HW 7
RK_4_A - HW 8
RK_4_B - HW 8

"""

"""
R_K_2_B C 
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

def RK4(a:float, b:float, N:int, x0:float, y0:float, actual_function=actual_function,
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


def RK2_B():
    """Compute the error ratios with taylor order 2 method"""
    h_list = [0.2, 0.1, 0.05, 0.025]
    N_list = [10, 20, 40, 80]
    y0 = 0.5
    x0 = 0.0
    b = 2.0
    a = 0.0

    x_list = []
    y_list = []
    #bounded_errors_list = []
    errors_list = []

    for h, N in zip(h_list, N_list):
        x, y, actual_vals = RK2(a, b, N, x0, y0)
        errors = taylor_high_order.compute_error(y, actual_vals)
        x_list.append(x)
        y_list.append(y)
        #bounded_errors_list.append(bounded_errors)
        errors_list.append(errors)

    #compute error ratios 
    r_list = []
    error_2 = taylor_high_order.get_every_nth_element(errors_list[1], 2)[1:]
    error_3 = taylor_high_order.get_every_nth_element(errors_list[2], 4)[1:]
    error_4 = taylor_high_order.get_every_nth_element(errors_list[3], 8)[1:]

    assert(error_2[-1] == errors_list[1][-1])
    assert(len(error_2) == len(error_3) == len(error_4))

    #divide errors
    np.seterr(divide='ignore', invalid='ignore')
    r1 = np.divide(errors_list[0][1:], error_2)
    r2 = np.divide(error_2, error_3)
    r3 = np.divide(error_3, error_4)

    #need to invert r1, r2, r3
    r1 = 1/r1
    r2 = 1/r2
    r3 = 1/r3
    
    return x, r1, r2, r3


def RK4_B():
    """Compute the error ratios with taylor order 2 method"""
    h_list = [0.2, 0.1, 0.05, 0.025]
    N_list = [10, 20, 40, 80]
    y0 = 0.5
    x0 = 0.0
    b = 2.0
    a = 0.0

    x_list = []
    y_list = []
    #bounded_errors_list = []
    errors_list = []

    for h, N in zip(h_list, N_list):
        x, y, actual_vals = RK4(a, b, N, x0, y0)
        errors = taylor_high_order.compute_error(y, actual_vals)
        x_list.append(x)
        y_list.append(y)
        #bounded_errors_list.append(bounded_errors)
        errors_list.append(errors)

    #compute error ratios 
    r_list = []
    error_2 = taylor_high_order.get_every_nth_element(errors_list[1], 2)[1:]
    error_3 = taylor_high_order.get_every_nth_element(errors_list[2], 4)[1:]
    error_4 = taylor_high_order.get_every_nth_element(errors_list[3], 8)[1:]

    assert(error_2[-1] == errors_list[1][-1])
    assert(len(error_2) == len(error_3) == len(error_4))

    #divide errors
    np.seterr(divide='ignore', invalid='ignore')
    r1 = np.divide(errors_list[0][1:], error_2)
    r2 = np.divide(error_2, error_3)
    r3 = np.divide(error_3, error_4)

    #need to invert r1, r2, r3
    # r1 = 1/r1
    # r2 = 1/r2
    # r3 = 1/r3
    
    return x, r1, r2, r3

if __name__ =="__main__":
    # R_K2_A()
    x, r1, r2, r3 = RK2_B()

    x4, r14, r24, r34 = RK4_B()



    