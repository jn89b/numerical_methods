import numpy as np
import pandas as pd

f = lambda x,y,yd : 4 + (0.25*x**3) - y * (yd/8)
fyv = lambda x,y,yd : -yd/8
fyvprime = lambda x,y,yd : -y/8


def gaussian_elimination(c,d,e,b):
    n = len(d)
    w = np.zeros(n)
    for i in range(1,n):
        s = -c[i]/d[i-1]
        d[i] = d[i] + s*e[i-1]
        b[i] = b[i] + s*b[i-1]
        # m = c[i-1]/d[i-1]
        # d[i] = d[i] - m*e[i-1]
        # b[i] = b[i] - m*b[i-1]

    W[n-1] = b[n-1]/d[n-1]
    
    for i in range(n-2,-1,-1):
        w[i] = (b[i] - e[i] * w[i+1])/d[i]
    
    return w


def F(x,w,alpha,beta,h):
    """F is the function that we want to find the root of."""
    n = len(w)
    fw = np.zeros(n)

    dy_1 = w[0]*((w[1]-alpha)/(16*h))
    fw[0] = -alpha +2*w[0] - w[1] + h**2 * (4 + (0.25 *x[0]**3) - dy_1)
    #fw[0] = -alpha +2*w[0] - w[1] + h**2 * f(x[0],w[0],(w[1]-alpha)/(2*h))

    for i in range(1,n-1):
        dy_i = (w[i]*(w[i+1]-w[i-1])/(16*h))
        fw[i] = -w[i-1] + 2*w[i] - w[i+1] + h**2 * (4 + (0.25*x[i]**3) - dy_i)
        # fw[i] = -w[i-1] + 2*w[i] - w[i+1] + h**2 * f(x[i],w[i],(w[i+1]-w[i-1])/(2*h))

    dy_n = w[n-1]*((beta - w[n-2])/(16*h))
    fw[n-1] = -w[n-2] + (2*w[n-1]) - beta + h**2 * (4 + (0.25*x[n-1]**3) - dy_n)
    # fw[n-1] = -w[n-2] + (2*w[n-1]) - beta + h**2 * f(x[n-1],w[n-1],(beta-w[n-2])/(2*h))

    return fw 

def compute_jacobian(w,alpha,beta,h):
    n = len(w)
    c = np.zeros(n)
    d = np.zeros(n)
    e = np.zeros(n)

    d[0] = 2 - (h/16) * ((w[1]-alpha))
    e[0] = -1 - h * ((w[0])/16)
    # d[0] = 2 + h**2 * fyv(x[0],w[0],(w[1]-alpha)/(2*h)) 
    # e[0] = -1 + h**2 * fyvprime(x[0],w[0],(w[1]-alpha)/(2*h)) * (1/(2*h))

    for i in range(1,n-1):
        c[i] = -1 + (h/16) * w[i]
        d[i] = 2 - h * ((w[i+1]-w[i-1])/16)
        e[i] = -1 - h * (w[i]/16)

        # d[i] = 2 + h**2 * fyv(x[i],w[i],(w[i+1]-w[i-1])/(2*h)) 
        # c[i] = -1 + h**2 * fyvprime(x[i],w[i],(w[i+1]-w[i-1])/(2*h)) * -(1/(2*h))
        # e[i] = -1 + h**2 * fyvprime(x[i],w[i],(w[i+1]-w[i-1])/(2*h)) * (1/(2*h))

    c[n-1] = -1 + h * (w[n-1]/16)
    d[n-1] = 2 - h * ((beta - w[n-2])/16)
    # d[n-1] = 2 + h**2 * fyv(x[n-1],w[n-1],(beta-w[n-2])/(2*h)) 
    # c[n-1] = -1 + h**2 * fyvprime(x[n-1],w[n-1],(beta-w[n-2])/(2*h)) * -(1/(2*h))
    
    return c,d,e


if __name__ == "__main__":
    a = 1 
    b = 3

    alpha = 17
    beta = 43/3
    n = 19
    h = (b-a)/(n+1)
    x = np.zeros(n)

    for i in range(n):
        x[i] = a + (i * h)+h
        
    W = np.zeros(n)
    #find linear line between a and b
    for i in range(n):
        W[i] = alpha + (alpha-beta)/(a-b) * (x[i] - a)

    i = 1
    Fw = F(x,W,alpha,beta,h)
    max_fw = max(abs(Fw))

    MF_array = []
    MF = max_fw
    MF_array.append(MF)

    tolerance = 10E-8
    max_iteration = 10

    while max_fw > tolerance:
        i += 1
        if i > max_iteration:
            print("i", i)
            break
        
        c,d,e = compute_jacobian(W,alpha,beta,h)
        P = gaussian_elimination(c,d,e,-Fw)
        W = W + P
        Fw = F(x,W,alpha,beta,h)
        max_fw = max(abs(Fw))
        MF = max_fw
        MF_array.append(MF)

    print(W)
    actual_w = [
    15.754503,
    14.771740,
    13.995677,
    13.386297,
    12.914252,
    12.557538,
    12.299326,
    12.126529,
    12.028814,
    11.997915,
    12.027142,
    12.111020,
    12.245025,
    12.425388,
    12.648944,
    12.913013,
    13.215312,
    13.553885,
    13.927046,
    14.333333
    ]
