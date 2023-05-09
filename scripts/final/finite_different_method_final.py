import numpy as np
import pandas as pd
import math 
p = lambda x : 2
q = lambda x : -1
r = lambda x : x*np.exp(x) - x


# c2 = -0.0392070132
# c1 = 1.1392070132
true_function = lambda x : x**3*np.exp(x)/6 - 5*x*np.exp(x)/3 + (np.exp(x)*2)- x - 2


def set_pandas_display_options() -> None:
    """Set pandas display options."""
    # Ref: https://stackoverflow.com/a/52432757/
    display = pd.options.display

    display.max_columns = 10000
    display.max_rows = 1000
    display.max_colwidth = 199
    display.width = 100

def gaussian_elimination(a,b,c,d,n):
    
    l = np.zeros((n+1))
    u = np.zeros((n+1))
    z = np.zeros((n+1))
    w = np.zeros((n+1))

    l[0] = a[0]
    u[0] = b[0]/l[0]
    z[0] = d[0]/l[0]

    for i in range(1,n):
        l[i] = a[i] - c[i] * u[i-1]
        u[i] = b[i]/l[i]
        z[i] = (d[i] - c[i] * z[i-1])/l[i]

    for j in range(n-1,-1,-1):
        w[j] = z[j] - (u[j] * w[j+1])
    return w
    

def finite_diff_method(aa,bb,alpha,beta,n, set_h=True):
    
    if set_h == True:
        h = (bb-aa)/(n+1)
    a = np.zeros((n+1))
    b = np.zeros((n+1))
    c = np.zeros((n+1))
    d = np.zeros((n+1))

    x_list = []
    x = aa + h
    x_list.append(x) 
    a[0] = 2 + h**2 * q(x)
    b[0] = -1 + 0.5 * h * p(x)
    d[0] = -h**2 * r(x) + (1 + 0.5 * h * p(x)) * alpha
    
    for i in range(1,n):
        x = x_list[-1] + h 
        x_list.append(x)
        a[i] = 2 + h**2 * q(x)
        c[i] = -1 - 0.5 * h * p(x)
        d[i] = -h**2 * r(x)
        
    for i in range(1,n-1):
        x = x_list[i]
        b[i] = -1 + 0.5 * h * p(x)
        
    x = bb - h
    x_list.append(x)
    # a[n] = 2 + h**2 * q(x)
    # c[n] = -1 - 0.5 * h * p(x)
    d[n-1] = -h**2 * r(x) + (1 - 0.5 * h * p(x)) * beta

    w = gaussian_elimination(a,b,c,d,n)
    return x_list,w

def problem19():
    aa = 0 
    bb = 2

    alpha = 0
    beta = -4
    n = 9

    x_list,w = finite_diff_method(aa,bb,alpha,beta,n)

    y_list = []
    error_list = []

    for estimated,val in zip(w,x_list):
        true_val = true_function(val)
        y_list.append(true_val)
        error_val = abs(true_val - estimated)
        error_list.append(error_val)

    # x_list = x.tolist()
    w_info = w.tolist()
    info_dict = {'x':x_list,
                 'y':y_list,
                 'w':w,
                 'error':error_list}

    df = pd.DataFrame(info_dict)
    
    print(df)

if __name__ == "__main__":
    #print in scientific notation for error
    pd.options.display.float_format = '{:,.4e}'.format

        
    set_pandas_display_options()    

    problem19()

    aa = 0 
    bb = 2

    alpha = 0
    beta = -4
    n = 9

    x_list_1,w_1 = finite_diff_method(aa,bb,alpha,beta,n)
    x_list_2,w_2 = finite_diff_method(aa,bb,alpha,beta,19)

    y_list_1 = []
    error_list_1 = []

    y_list_2 = []
    error_list_2 = []

    for estimated,val in zip(w_1,x_list_1):
        true_val = true_function(val)
        y_list_1.append(true_val)
        error_val = abs(true_val - estimated)
        error_list_1.append(error_val)

    for estimated,val in zip(w_2,x_list_2):
        true_val = true_function(val)
        y_list_2.append(true_val)
        error_val = abs(true_val - estimated)
        error_list_2.append(error_val)

    #w_1 to array
    error1 = np.array(error_list_1)
    error2 = np.array(error_list_2)
    error2 = error2[1::2]

    ratio = error1/error2
    # print(ratio[:-1])

    last_w1 = w_1[-2]
    last_w2 = w_2[-2]

    w3 = (4*w_2[1::2] - w_1)/3
    w_2 = w_2[1::2]

    estimated_error1_list = []
    estimated_error2_list = []

    for estimated1, estimated2, approx_true in zip(w_1,w_2,w3):
        estimated_error1 = abs(approx_true - estimated1)
        estimated_error2 = abs(approx_true - estimated2)
        estimated_error1_list.append(estimated_error1)
        estimated_error2_list.append(estimated_error2)

    info_dict = {'x':x_list_1[:-1],
                 'true error': error1[:-1],
                 'E1E':estimated_error1_list[:-1]}
    
    max_estimated_error = max(estimated_error1_list[:-1])
    c = max_estimated_error/(0.2)**2

    desired_error = 10E-4
    h_desired = np.sqrt(desired_error/c)
    print(h_desired)
    
    #round h_desired to nearest 0.1
    
    df = pd.DataFrame(info_dict)

    #print in scientific notation for error

    print(df)
    n_desired = math.ceil(((bb-aa)/h_desired) - 1)    
    x_final, w_final  = finite_diff_method(aa,bb,alpha,beta,n_desired)
    
    y_list = []
    error_list = []
    
    for estimated,val in zip(w_final,x_final):
        true_val = true_function(val)
        y_list.append(true_val)
        error_val = abs(true_val - estimated)
        error_list.append(error_val)

    print("c = ",c)
    print("h_desired = ",h_desired)
    print("n_desired = ",n_desired)
    print("max estimated error = ", max(error_list[:-1]))

    # x_list = x.tolist()
    w_info = w_final.tolist()
    info_dict = {'x':x_final,
                 'y':y_list,
                 'w':w_final,
                 'error':error_list}

    df = pd.DataFrame(info_dict)
    
    print(df)




