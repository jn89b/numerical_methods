import numpy as np
import pandas as pd

p = lambda x : -2/x
q = lambda x : (2 / x**2)
r = lambda x : (np.sin(np.log(x)))/(x**2)


c2 = -0.0392070132
c1 = 1.1392070132
# true_function = lambda x : (c1*x) + (c2/(x**2)) \
#     - ((0.3)*np.sin(np.log(x))) \
#     - ((0.1)*np.cos(np.log(x)))  

true_function = lambda x : x**2 + (16/x)

f = lambda x,y,y_prime: 4 + (0.25*x**3) - ((y*y_prime)/8)
fy = lambda x,y,y_prime: y_prime/8
fy_prime = lambda x,y,y_prime: - y/8


def set_pandas_display_options() -> None:
    """Set pandas display options."""
    # Ref: https://stackoverflow.com/a/52432757/
    display = pd.options.display

    display.max_columns = 10000
    display.max_rows = 10000
    display.max_colwidth = 199
    display.width = 1000

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
    a = np.zeros((n+1))
    b = np.zeros((n+1))
    c = np.zeros((n+1))
    d = np.zeros((n+1))
    w = np.zeros((n+1))

    h = (bb-aa)/(n)

    x_list = []
    x = aa + h
    x_list.append(x) 

    w[0] = alpha
    w[-1] = beta
    for i in range(n):
        w[i] = aa + (i*h*((bb-aa)/(bb-x_list[-1])))
        x = x_list[-1] + h 
        x_list.append(x)

    t = (w[1] - alpha)/ (2*h)
    a[0] = 2 + h**2 * fy(x, w[0], t)
    b[0] = -1 + (h/2) * fy_prime(x, w[0], t)
    d[0] = -((2*w[0]) -w[1] -aa + (h**2 * f(x, w[0], t)))

    for i in range(1,n):
        t = (w[i+1] - w[i-1])/(2*h)
        x = x_list[i]
        a[i] = 2 + h**2 * fy(x, w[i], t)
        c[i] = -1 - (h/2) * fy_prime(x, w[i], t)
        d[i] = -((2*w[i]) -w[i+1] -w[i-1] -aa + (h**2 * f(x, w[i], t)))

    for i in range(1,n-1):
        x = x_list[i]
        b[i] = -1 + (h/2) * fy_prime(x, w[i], t)
        
    x = bb - h
    x_list.append(x)
    d[n-1] = -((2*w[i]) -w[i+1] -w[i-1] -aa + (h**2 * f(x, w[i], t)))

    w = gaussian_elimination(a,b,c,d,n)

    return x_list,w

if __name__ == "__main__":
    aa = 1 
    bb = 3
    a = 17
    b = 43/35
    n = 20
    x_list,w = finite_diff_method(aa,bb,alpha,beta,n)

    true_y_list = []
    for x in x_list:
        true_y = true_function(x)   
        true_y_list.append(true_y)

    dict_info = {
        "x_list": x_list,
        "true_y_list": true_y_list
    }

    df = pd.DataFrame(dict_info)
    set_pandas_display_options()