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

another_function = lambda x, y: (2*x*y)/(x**2 - y**2)

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

def RK4(a:float, b:float, N:int, x0:float, y0:float, 
    actual_function=actual_function,
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

        #print in high precision
        print("s1 = ", s1, "s2 = ", s2, "s3 = ", s3, "s4 = ", s4)
        
        y[i+1] = y[i] + (h*(s1 + 2*s2 + 2*s3 + s4)/6)
        print(y[i+1])
        actual_vals[i+1] = actual_function(x[i+1])

    return x, y, actual_vals

def R_K2_A():
    """RK2A Due Next Wednesday """
    a = 0 
    b = 2
    ad = 0.5
    N = 10

    x, y, actual_vals = RK2(a, b, N, a, ad)
    error = y - actual_vals
    print("RK2A")
    print("x\t y\t actual\t error")

    for i in range(N+1):

        #keep 4 decimal places
        x[i] = round(x[i], 2)
        y[i] = round(y[i], 4)
        actual_vals[i] = round(actual_vals[i], 4)
        error[i] = round(error[i], 4)

        print(x[i], "\t", y[i], "\t", actual_vals[i], "\t", error[i])



def RK4_A():
    """RK4A Due Next Wednesday """
    a = 0 
    b = 2
    ad = 0.5
    N = 10

    x, y, actual_vals = RK4(a, b, N, a, ad)
    error = y - actual_vals

    # #print in table format
    # for i in range(N+1):
    #     print("x = ", x[i], "y = ", y[i], "actual = ", actual_vals[i])
    print("RK4A")
    print("x\t y\t actual\t error")

    for i in range(N+1):

        #keep 4 decimal places
        x[i] = round(x[i], 2)
        y[i] = round(y[i], 4)
        actual_vals[i] = round(actual_vals[i], 4)

        #error scientific notation
        error[i] = round(error[i], 8)
        #error[i] = round(error[i], 4)

        print(x[i], "\t", y[i], "\t", actual_vals[i], "\t", error[i])


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
    errors_list = []

    for h, N in zip(h_list, N_list):
        x, y, actual_vals = RK2(a, b, N, x0, y0)
        errors = taylor_high_order.compute_error(y, actual_vals)
        x_list.append(x)
        y_list.append(y)
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
    print("x = ", x_list[0])
    #print in table format
    x_list = x_list[0][1:]
    y_list = y_list[0][1:]
    # errors_list = errors_list[:][1:]
    
    print("h\t x\t y\t r1\t r2\t r3")
    for i in range(len(r1)):
        h_list[0] = round(h_list[0], 2)
        
        #round x_list to 2 decimal places
        x_list[i] = round(x_list[i], 2)
        y_list[i] = round(y_list[i], 4)
        # y_list[0][i] = round(y_list[i], 4)
        actual_vals[i] = round(actual_vals[i], 4)
        errors_list[0][i] = round(errors_list[0][i], 4)
        r1[i] = round(r1[i], 4)
        r2[i] = round(r2[i], 4)
        r3[i] = round(r3[i], 4)

        print(h_list[0], "\t", x_list[i], "\t", y_list[i], 
        "\t", r1[i], "\t", r2[i], "\t", r3[i])

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

    x_list = x_list[0][1:]
    y_list = y_list[0][1:]
    # errors_list = errors_list[:][1:]
    
    print("h\t x\t y\t r1\t r2\t r3")
    for i in range(len(r1)):
        h_list[0] = round(h_list[0], 2)
        
        #round x_list to 2 decimal places
        x_list[i] = round(x_list[i], 2)
        y_list[i] = round(y_list[i], 4)
        # y_list[0][i] = round(y_list[i], 4)
        actual_vals[i] = round(actual_vals[i], 4)
        errors_list[0][i] = round(errors_list[0][i], 4)
        r1[i] = round(r1[i], 4)
        r2[i] = round(r2[i], 4)
        r3[i] = round(r3[i], 4)

        print(h_list[0], "\t", x_list[i], "\t", y_list[i], 
        "\t", r1[i], "\t", r2[i], "\t", r3[i])

    return x, r1, r2, r3


def Problem10A():
    """Compute the RK4 method"""
    x1,w1, y1 = RK4(a=0, b=2, N=20, x0=0, y0=3, function=another_function)
    x2,w2, y2 = RK4(a=0, b=2, N=40, x0=0, y0=3, function=another_function)
    x3,w3, y3 = RK4(a=0, b=2, N=80, x0=0, y0=3, function=another_function)

    print("w1", w1[-1])
    print("w2", w2[-1])
    print("w3", w3[-1])

if __name__ =="__main__":
    # R_K2_A()
    # x, r1, r2, r3 = RK2_B()

    # RK4_A()
    # x4, r14, r24, r34 = RK4_B()

    a = 0 
    b = 0.2
    y0 = 0.5
    N = 1
    
    Problem10A()
        

    