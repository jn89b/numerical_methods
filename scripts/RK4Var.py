"""
Runge Kutta Variable Step Size

12A - 27 in 5.4

12B 28 in 5.4: 
scatter plot estimate, actual along time

Outputs should be:
    Step Size has become too small at N = 36

Format as long:
0.1         0.1     7.99978     
0.3         0.2     7.99936     
...
1505.9      0.4     0.24..  
1506.1      0.2     0.1204..
1506.125    0.4
1506.1375   
...
1506.141796875  0.000390625  0.008891

HW13A:
function = 2*t*y / (t**2 - y**2)
a = 0;
b = 2;
al = 3;
tol = 10E-4;
[T,H,W] = RK4Var(a,b,al,tol);
N = length(T);

plot(T,W,'o',T,W,'-');
Outputs should be:
    t = 2 -
    t = 2 -
    t = 2 - way different than previous two approximations (destabilization)
    Step-size has become too small at N = 16

    THW:
    0.1     0.1            
    0.3     0.2
    0.7     0.4
    ...
    1.4996  0.00039063
    1.4998  0.00019351 

PLOT this data to visualize what will happen

HW13B:
Verify that :
    t^2 + y^2 = ky

Satisfies the DE, y' = f(t,y) = 2*t*y / (t^2 - y^2) for any k

Determine k in t^2 + y^2 = ky
so that y(0) = 3, k=3 is the solution to the DE

Now determine the largest possible b (end time) 
for on which t*2 + y^2 = 3y is a solution to the DE, 
before the solution becomes unstable

Numerical solution:
    t^2 + y^2 = 3y is a circle so we can complete the square to get:
    t^2 + y^2 -3y = 0
    (t)^2 + (y-1.5)^2 = 1.5^2

    derivative of this is:
    So derivative of the value is not defined 
    
"""

import numpy as np 
import matplotlib.pyplot as plt
import pandas as pd

r = 0.1
g = 32.1 
f = lambda x,y : -0.6*np.pi*r**2*np.sqrt(2*g)*np.sqrt(y)/ (np.pi*y**2)

k = 6.22E-19
n_1 = 2E3
n_2 = 2E3
n_3 = 3E3 
chemical_function = lambda x,y:\
    k*(n_1-y/2)**2*(n_2-(y/2))**2*(n_3-(3*y/4))**3

thirteen_a_function = lambda x,y: 2*x*y / (x**2 - y**2)

true_current1 = lambda x: -3.375*np.exp(-2*x) + \
    1.875*np.exp(-0.4*x) + 1.5

true_current2 = lambda x: -2.25*np.exp(-2*x) + \
    2.25*np.exp(-0.4*x)


true_val1 = lambda x: 0.2*np.exp(2*x) * (np.sin(x) - (2*np.cos(x)))
true_val2 = lambda x: 0.2*np.exp(2*x) * (4*np.sin(x) - (3*np.cos(x)))

def set_pandas_display_options() -> None:
    """Set pandas display options."""
    # Ref: https://stackoverflow.com/a/52432757/
    display = pd.options.display

    display.max_columns = 1000
    display.max_rows = 1000
    display.max_colwidth = 199
    display.width = 1000
    # display.precision = 2  # set as needed

def RK4(N:int, x0:float, y0:float, 
        h = 0.1,
        function=f):
    """Computes RK4 for a given function"""

    w = np.zeros(N+1)
    x = np.zeros(N+1)
    x[0] = x0
    w[0] = y0

    for i in range(N):
        x[i+1] = x[i] + h
        s1 = function(x[i], w[i])
        s2 = function(x[i] + h/2, w[i] + h*s1/2)
        s3 = function(x[i] + h/2, w[i] + h*s2/2)
        s4 = function(x[i] + h, w[i] + h*s3)

        w[i+1] = w[i] + (h*(s1 + 2*s2 + 2*s3 + s4)/6)

    return x,w

def AdaptiveRK(a:float, b:float, y0:float,
               h_init = 0.1, 
               min_h = 10E-4,
               tolerance = 0.1,
               function = chemical_function):
    
    T = []
    W = []
    H = []
    t = a 
    w = y0
    h = h_init 
    
    if a + h > b:
        h = b - a
    
    bound_tolerance = b - (tolerance)

    while t < bound_tolerance:

        x1,w1 = RK4(1, t, w, h, function=function)
        x2, w2 = RK4(2, t, w, h/2, function=function)
        w3 = (16*w2[-1] - w1[-1])/15

        error = abs(w3 - w1[-1])

        if error < tolerance:
            w = w3
            t = t + h
            T.append(t)
            W.append(w)
            H.append(h)

            if error < tolerance/128:
                h = 2*h
            
            if t + h > b:
                h = b - t
        else:
            h = h/2
            if h < min_h:
                print("Step size has become too small at N = ", len(T))
                return T,W,H
        
    return T,W,H

def Problem12B():
    a = 0.0#0
    b = 1600
    al = 8
    tol = 10E-5
    min_h = 10E-5

    T,W,H = AdaptiveRK(a,b,al,
                       tolerance=tol, 
                       min_h=min_h, 
                       function=f)

    info_dict = {'t':T, 'h':H, 'w':W}
    df = pd.DataFrame(info_dict)
    print(df)
    N = len(T)

    plt.plot(T,W,'o',T,W,'-')
    plt.show()

def Problem13A():
    a = 0
    b = 2
    al = 3
    tol = 10E-5
    min_h = 10E-5
    [T,W,H] = AdaptiveRK(a,b,al,
                         tolerance=tol,
                         min_h=min_h,
                         function=thirteen_a_function)

    N = len(T)
    df = pd.DataFrame({'t':T, 'h':H, 'w':W})
    print(df)
    plt.plot(T,W,'o',T,W,'-')
    plt.show()
    

def RK4Vector(a,b, w, N, h=[False,0.1], condition= 1):
    """
    RK4 for a vector function
    """

    T = np.zeros(N+1)
    W = np.zeros((N+1, len(w)))
    
    W[0] = w
    T[0] = a
    
    if h[0] == False:
        h[1] = (b-a)/N
    else:
        h[1] = h[1]

    h_val = h[1]

    if condition == 1:
        function = vector_function
    else:
        function = higher_order_function

    for i in range(N):
        s1 = function(T[i], W[i])
        s2 = function(T[i] + h_val/2, W[i] + h_val*s1/2)
        s3 = function(T[i] + h_val/2, W[i] + h_val*s2/2)
        s4 = function(T[i] + h_val, W[i] + h_val*s3)
        W[i+1] = W[i] + (h_val*(s1 + 2*s2 + 2*s3 + s4)/6)
        T[i+1] = T[i] + h_val
    
    return T,W
        
def Problem14A():
    a = 0
    b = 10
    w = np.array([0,0])
    N = 100
    
    T,W = RK4Vector(a,b,w,N)
    W1 = W[:,0]
    W2 = W[:,1]
    y1 = true_current1(T)
    y2 = true_current2(T)
    E1 = abs(y1-W1)
    E2 = abs(y2-W2)

    info_dict = {
        't':T, 
        'w1':W1, 
        'w2':W2, 
        'E1':E1, 
        'E2':E2
    }
    
    df = pd.DataFrame(info_dict)
    print(df)
    

def vector_function(t,w):
    f1 = lambda x,w: (-4*w[0]) + (3*w[1]) + 6
    f2 = lambda x,w: (-2.4*w[0]) + (1.6*w[1]) + 3.6 
    
    fv = np.zeros(len(w))
    fv[0] = f1(t,w)
    fv[1] = f2(t,w)
    return fv

def AdaptiveRKVector(a,b,w,tol,condition=1):
    T = []
    W = []
    H = []
    t = a 
    w = w
    h = 0.1
    
    if a + h > b:
        h = b - a
    
    while t < b - 10E-12:
        x1,w1 = RK4Vector(t, b, w,1, h=[True, h], condition=condition)#RK4Vector(t, w, h, N)
        x2, w2 = RK4Vector(t, b, w, 2, h=[True,h/2], condition=condition)
        last_w1 = w1[-1]
        last_w2 = w2[-1]

        w3 = (16*last_w2 - last_w1)/15  
        error = max(abs(w3 - last_w1))

        if error < tol:
            w = w3
            t = t + h
            T.append(t)
            W.append(w)
            H.append(h)

            if error < tol/128:
                h = 2*h
            
            if t + h > b:
                h = b - t

        else:
            h = h/2
            if h < 10E-4:
                print("Step size has become too small at N = ", len(T))
                return T,W,H

    return T,W,H

def Problem14B():
    a = 0
    b = 10
    w = np.array([0,0])
    N = 100

    T,W,H = AdaptiveRKVector(a,b,w,10E-5)
    T = np.array(T)
    W = np.array(W)
    N = len(T)
    W1 = W[:,0]
    W2 = W[:,1]
    y1 = true_current1(T)
    y2 = true_current2(T)
    E1 = abs(y1-W1)
    E2 = abs(y2-W2)

    info_dict = {
        't':T,
        'w1':W1,
        'w2':W2,
        'E1':E1,
        'E2':E2
    }

    df = pd.DataFrame(info_dict)
    print(df)


def higher_order_function(t,w):
    f1 = lambda x,w: w[1]
    f2 = lambda x,w: -2*w[0] + 2*w[1] + (np.exp(2*x)*np.sin(x))
    
    fv = np.zeros(len(w))
    fv[0] = f1(t,w)
    fv[1] = f2(t,w)
    return fv

def Problem15A():
    a = 0 
    b = 3
    w = np.array([-0.4,-0.6])
    tol = 10E-5
    N = 30
    #T,W,H = AdaptiveRKVector(a,b,w,tol, condition=2)
    T,W = RK4Vector(a,b,w,N, condition=2)
    T = np.array(T)
    W = np.array(W)
    #N = len(T)
    W1 = W[:,0]
    W2 = W[:,1]
    y1 = true_val1(T)
    y2 = true_val2(T)
    E1 = abs(y1-W1)
    E2 = abs(y2-W2)

    info_dict = {
        't':T,
        'w1':W1,
        'w2':W2,
        'E1':E1,
        'E2':E2
    }

    df = pd.DataFrame(info_dict)
    print(df)    

def Problem15B():
    a = 0 
    b = 3
    w = np.array([-0.4,-0.6])
    tol = 10E-5
    T,W,H = AdaptiveRKVector(a,b,w,tol, condition=2)
    T = np.array(T)
    W = np.array(W)
    N = len(T)
    W1 = W[:,0]
    W2 = W[:,1]
    y1 = true_val1(T)
    y2 = true_val2(T)
    E1 = abs(y1-W1)
    E2 = abs(y2-W2)

    info_dict = {
        't':T,
        'w1':W1,
        'w2':W2,
        'E1':E1,
        'E2':E2
    }

    df = pd.DataFrame(info_dict)
    print(df)


if __name__ == "__main__":
    set_pandas_display_options()
    #Problem12B()
    #Problem13A()
    
    #Problem14A()
    #Problem14B()
    
    #Problem15A()
    #Problem15B()

    #Problem
    

        
    



