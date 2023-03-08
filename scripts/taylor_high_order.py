import numpy as np
import math as m
import decimal

actual_function = lambda x: (x+1)*(x+1) - 0.5*m.exp(x)
approx_function = lambda x, y: y - x**2 + 1
approx_dy = lambda x, y: approx_function(x, y) - 2*x
approx_ddy = lambda x, y: approx_function(x, y) - 2*x - 2
approx_dddy = lambda x, y: approx_function(x, y) - 2*x - 2
approx_ddddy = lambda x, y: approx_function(x, y) - 2*x - 2

EULER_NUMBER = 2.71828182845904523536


def compute_bound_error_estimate(h:float, ME:float,
    x:float, factor_number:float, raise_number:float):

    return ME* (EULER_NUMBER**(x) - 1) * (h**raise_number) / factor_number 

def compute_best_h(error_tol:float, ME:float, 
    x:float, factor_number:float, raise_number:float):

    return ((error_tol * factor_number) / (ME * (EULER_NUMBER**(x) - 1)))**(1/raise_number)



def taylor_method4(a:float, b:float, N:int, 
    x0:float, y0:float, h=None, compute_bound_error=False,
    find_best_h=False, order=4):
    
    """Computes Taylor Method for a given function"""
    x = np.zeros(N+1)
    y = np.zeros(N+1)
    actual_vals = np.zeros(N+1)

    if compute_bound_error:
        error_estimated_margins = np.zeros(N+1)
        bounded_estimated_error = []

    x[0] = x0
    y[0] = y0
    actual_vals[0] = actual_function(x0)
    
    #round to ceiling of h
    if h is None:
        h = (b-a)/N

    for i in range(N):
        x[i+1] = x[i] + h

        if order==4:
            y[i+1] = y[i] + (h*approx_function(x[i], y[i])) \
                + (h**2*approx_dy(x[i], y[i]))/2 \
                    + (h**3*approx_ddy(x[i], y[i]))/6 \
                        + (h**4*approx_dddy(x[i], y[i]))/24

        if order==2:
            y[i+1] = y[i] + (h*approx_function(x[i], y[i])) \
                + (h**2*approx_dy(x[i], y[i]))/2

        if compute_bound_error:
            error_estimated_margins[i+1] = abs(approx_ddddy(x[i+1], y[i+1]))

        actual_vals[i+1] = actual_function(x[i+1])

    if order == 4:
        raise_number = 4
        factor_number = 120

    if order == 2:
        raise_number = 2
        factor_number = 6

    if compute_bound_error:
        max_estimated_error = max(error_estimated_margins)

        for xi in x:
            bounded_estimated_error.append(abs(compute_bound_error_estimate(
                h, max_estimated_error, xi, factor_number, raise_number)))
            
    if find_best_h:
        ME = max_estimated_error
        best_h =  compute_best_h(1e-4, ME, x[-1], factor_number, raise_number)

        return x,y, actual_vals, bounded_estimated_error, best_h

    if compute_bound_error:
        return x, y, actual_vals, bounded_estimated_error

    return x, y, actual_vals
 

def compute_error(y, actual_vals):
    """"""
    return np.abs(y - actual_vals)

def get_every_nth_element(array, n):
    #remove the first element
    return array[::n]

def problem_5_a_b():

    h = 0.2 
    y0 = 0.5 
    x0 = 0.0
    N = 10 
    b = 2.0
    a = 0.0

    x, y, actual_vals, bounded_errors = taylor_method4(a, b, N, x0, y0, 
        compute_bound_error=True)
    errors = compute_error(y, actual_vals)

    return x, y, actual_vals, errors, bounded_errors


def problem_4_a():
    """Compute the error ratios with taylor order 2 method"""
    h_list = [0.2, 0.1, 0.05, 0.025]
    N_list = [10, 20, 40, 80]
    y0 = 0.5
    x0 = 0.0
    b = 2.0
    a = 0.0

    x_list = []
    y_list = []
    bounded_errors_list = []
    errors_list = []

    for h, N in zip(h_list, N_list):
        x, y, actual_vals, bounded_errors = taylor_method4(a, b, N, x0, y0, 
            compute_bound_error=True, order=2)
        errors = compute_error(y, actual_vals)
        x_list.append(x)
        y_list.append(y)
        bounded_errors_list.append(bounded_errors)
        errors_list.append(errors)

    #compute error ratios 
    r_list = []
    error_2 = get_every_nth_element(errors_list[1], 2)[1:]
    error_3 = get_every_nth_element(errors_list[2], 4)[1:]
    error_4 = get_every_nth_element(errors_list[3], 8)[1:]

    assert(error_2[-1] == errors_list[1][-1])
    assert(len(error_2) == len(error_3) == len(error_4))

    #divide errors
    np.seterr(divide='ignore', invalid='ignore')
    r1 = np.divide(errors_list[0][1:], error_2)
    r2 = np.divide(error_2, error_3)
    r3 = np.divide(error_3, error_4)

    return x, r1, r2, r3

def problem_5_b():
    x,y, actual_vals, bounded_estimated_error, best_h = taylor_method4(
        a, b, N, x0, y0, h=h, compute_bound_error=True, find_best_h=True, 
        order=2)

    #round up to the nearest integer

    best_n = m.ceil((b-a)/best_h)
    x_best, y_best, actual_vals_best, bounded_estimated_error_best = taylor_method4(
        a, b, best_n, x0, y0, h=best_h, compute_bound_error=True,
        order=2)

    error_list = []
    for w0, y_actual in zip(y_best, actual_vals_best):
        error_list.append(abs(w0 - y_actual))
    
    for error in error_list:
        print(error)
    
    best_error = max(error_list)
    print(min(compute_error(y_best, actual_vals_best)))
    return best_n, best_error, best_h


def problem_6_a():
    h_list = [0.2, 0.1, 0.05, 0.025]
    N_list = [10, 20, 40, 80]
    y0 = 0.5
    x0 = 0.0
    b = 2.0
    a = 0.0

    x_list = []
    y_list = []
    bounded_errors_list = []
    errors_list = []

    for h, N in zip(h_list, N_list):
        x, y, actual_vals, bounded_errors = taylor_method4(a, b, N, x0, y0, 
            compute_bound_error=True)
        errors = compute_error(y, actual_vals)
        x_list.append(x)
        y_list.append(y)
        bounded_errors_list.append(bounded_errors)
        errors_list.append(errors)

    #compute error ratios 
    r_list = []
    error_2 = get_every_nth_element(errors_list[1], 2)[1:]
    error_3 = get_every_nth_element(errors_list[2], 4)[1:]
    error_4 = get_every_nth_element(errors_list[3], 8)[1:]

    assert(error_2[-1] == errors_list[1][-1])
    assert(len(error_2) == len(error_3) == len(error_4))

    #divide errors
    np.seterr(divide='ignore', invalid='ignore')
    r1 = np.divide(errors_list[0][1:], error_2)
    r2 = np.divide(error_2, error_3)
    r3 = np.divide(error_3, error_4)

    return x, r1, r2, r3


def problem_6_b():
    x,y, actual_vals, bounded_estimated_error, best_h = taylor_method4(
        a, b, N, x0, y0, h=h, compute_bound_error=True, find_best_h=True)

    #round up to the nearest integer
    best_n = m.ceil((b-a)/best_h)
    x_best, y_best, actual_vals_best, bounded_estimated_error_best = taylor_method4(
        a, b, best_n, x0, y0, h=best_h, compute_bound_error=True)

    best_error = max(compute_error(y_best, actual_vals_best))
    
    return best_n, best_error, best_h

if __name__ == "__main__":

    y0 = 0.5
    x0 = 0.0
    b = 2.0
    a = 0.0
    N = 10
    h = 0.2

    x,r1, r2, r3 = problem_4_a()
    x = get_every_nth_element(x, 8)
    x = [round(i, 5) for i in x]
    r1 = [round(i, 5) for i in r1]
    r2 = [round(i, 5) for i in r2]
    r3 = [round(i, 5) for i in r3]

    # #output x, r1, r2, r3 as table to text file
    # with open("problem_4_a.txt", "w") as f:
    #     f.write("x\t r1\t r2\t r3")
    #     for i in range(len(x)+1):

    #         if i == len(x):
    #             f.write("\n")
    #             f.write("x = 2.0\t r1 = 1.0\t r2 = 1.0\t r3 = 1.0")
    #             break

    #         f.write("\n")
    #         f.write(str(x[i]) + "\t" + str(r1[i]) + "\t" + \
    #             str(r2[i]) + "\t" + str(r3[i]))



    best_taylor2_n, best_taylor2_error, best_taylor2_h = problem_5_b()
    print("best_taylor2_n: ", best_taylor2_n)
    print("best_taylor2_error: ", best_taylor2_error)
    print("best_taylor2_h: ", best_taylor2_h)

    ## 5_a_b_c
    x,y,actual_vals,errors,bounded_errors = problem_5_a_b()
    x = [round(i, 4) for i in x]
    y = [round(i, 4) for i in y]
    actual_vals = [round(i, 4) for i in actual_vals]

    #display errors in scientific notation
    errors = ["{:.2e}".format(i) for i in errors]
    bounded_errors = ["{:.2e}".format(i) for i in bounded_errors]


    # #output x, y, actual_vals, errors, bounded_errors as table to text file
    with open("problem_5_a_b_c.txt", "w") as f:
        f.write("x\t y\t actual_vals\t errors\t bounded_errors")
        for i in range(len(x)+1):
                
                if i == len(x):
                    f.write("\n")
                    f.write("x = 2.0\t y = 0.5\t actual_vals = 0.5\t errors = 0.0\t bounded_errors = 0.0")
                    break
        
                f.write("\n")
                f.write(str(x[i]) + "\t" + str(y[i]) + "\t" + \
                    str(actual_vals[i]) + "\t" + str(errors[i]) + "\t" + str(bounded_errors[i]))


    test,r1, r2, r3 = problem_6_a()
    decimal_precision = 4

    x = [round(i, decimal_precision) for i in x]
    r1 = [round(i, decimal_precision) for i in r1]
    r2 = [round(i, decimal_precision) for i in r2]
    r3 = [round(i, decimal_precision) for i in r3]

    # #output x, r1, r2, r3 as table to text file
    # with open("problem_6_a.txt", "w") as f:
    #     f.write("x\t r1\t r2\t r3")
    #     for i in range(len(x)):
            
    #         if i == len(x):
    #             f.write("\n")
    #             f.write("x = 2.0\t r1 = 1.0\t r2 = 1.0\t r3 = 1.0")
    #             break
    #         else:
    #             f.write("\n")
    #             f.write(str(x[i]) + "\t" + str(r1[i]) + "\t" + \
    #                 str(r2[i]) + "\t" + str(r3[i]))



    best_n, best_error, best_h = problem_6_b()

    print("best_n: ", best_n)
    print("best_error: ", best_error)
    print("best_h: ", best_h)


