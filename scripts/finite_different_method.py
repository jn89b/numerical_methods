import numpy as np
import pandas as pd

p = lambda x : -2/x
q = lambda x : (2 / x**2)
r = lambda x : (np.sin(np.log(x)))/(x**2)


c2 = -0.0392070132
c1 = 1.1392070132
true_function = lambda x : (c1*x) + (c2/(x**2)) \
    - ((0.3)*np.sin(np.log(x))) \
    - ((0.1)*np.cos(np.log(x)))  


def set_pandas_display_options() -> None:
    """Set pandas display options."""
    # Ref: https://stackoverflow.com/a/52432757/
    display = pd.options.display

    display.max_columns = 1000
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
    

def finite_diff_method(aa,bb,alpha,beta,n):
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
    aa = 1 
    bb = 2

    alpha = 1
    beta = 2
    n = 9

    x_list,w = finite_diff_method(aa,bb,alpha,beta,n)

    y_list = []
    error_list = []

    for estimated,val in zip(w,x_list):
        true_val = true_function(val)
        y_list.append(true_val)
        error_val = abs(true_val - estimated)
        error_list.append(error_val)
        #print(true_function(val))f

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
    pd.options.display.float_format = '{:,.2e}'.format
        
    set_pandas_display_options()    

    problem19()

    aa = 1 
    bb = 2

    alpha = 1
    beta = 2
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
    
    df = pd.DataFrame(info_dict)
    #print in scientific notation for error
    pd.options.display.float_format = '{:,.5e}'.format

    print(df)

    



