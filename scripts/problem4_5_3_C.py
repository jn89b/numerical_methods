import numpy as np
import math 
import pandas as pd

fv = lambda x,y: y - x**2 + 1
fvd = lambda x,y: fv(x,y) - 2*x
fv2d = lambda x,y: fv(x,y) - 2*x - 2
fv3d = lambda x,y: fv(x,y) - 2*x - 2

actual_function = lambda x,y: (x+1)**2  - 0.5*math.exp(x)

def Taylor4(a:float, b:float, y0:float, N:int):
    
    h = (b-a)/N
    x = np.zeros(N+1)
    w = np.zeros(N+1)
    error = np.zeros(N+1)
    
    x[0] = a 
    w[0] = y0

    for i in range(N):

        x[i+1] = x[i] + h

        w[i+1] = w[i] + (h*fv(x[i], w[i])) \
            + (h**2*fvd(x[i], w[i]))/2 \
            + (h**3*fv2d(x[i], w[i]))/6 \
            + (h**4*fv3d(x[i], w[i]))/24
        
        error[i+1] = abs(actual_function(x[i+1], w[i+1]) - w[i+1])
        
    return h,x,w, error

def problem4_5_3_C():
    a = 0 
    b = 2 
    y0 = 0.5
    N = 10

    h, x, w, error= Taylor4(a, b, y0, N)    
    h1, x1, w1, error1= Taylor4(a, b, y0, 2*N)
    h2, x2, w2, error2 = Taylor4(a, b, y0, 4*N)
    h3, x3, w3, error3 = Taylor4(a, b, y0, 8*N)

    error1 = error1[::2]
    error2 = error2[::4]
    error3 = error3[::8]

    r1 = error/error1
    r2 = error1/error2
    r3 = error2/error3

    #create dataframe for table
    df = pd.DataFrame({'x':x, 'w':w, 'error':error, 'r1':r1, 'r2':r2, 'r3':r3})
    df = df.round(4)

    return df

df = problem4_5_3_C()
print(df)




