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

v = 0.1
g = 32.1
# some_function = lambda x,y: x**2 + y**2

some_function = lambda x,y:\
     2*x*y / (x**2 - y**2)
# some_function = lambda r,y: \
#     (-0.6*np.pi*r**2)*(np.sqrt(2*g)) * ((np.sqrt(y))/ (np.pi*y**2))

k = 6.22E-19
n_1 = 2E3
n_2 = 2E3
n_3 = 3E3 
chemical_function = lambda x,y:\
    k*(n_1-y/2)**2*(n_2-(y/2))**2*(n_3-(3*y/4))**3


def RK(t,w,h,n, some_function=chemical_function):
    for i in range(n):
        s1 = some_function(t, w)
        s2 = some_function(t + h/2, w + (h*s1)/2)
        s3 = some_function(t + h/2, w + (h*s2)/2)
        s4 = some_function(t + h, w + h*s3)
        print(s1,2*s2,2*s3,s4)
        w= w + ((s1 + 2*s2 + 2*s3 + s4)/6)
        t = t + h

    return t, w

def RK4Var(a, b, w, tol,
           h_init = 0.1,
           min_h = 10E-4, 
           some_function=chemical_function):
    T = []
    W = []
    H = []
    t = a
    h = h_init
    min_h = min_h

    if a + h > b:
        h = b - a

    while t < (b-10E-12):

        w1 = RK(t, w, h, 1, some_function)  
        w2 = RK(t, w, h/2, 2, some_function)
        w3 = (16*w2[0]- w1[0])/15
        print("w1", w1)
        print("w2", w2)
        error = abs(w3 - w1[0])
        print("w3: ", w3)
        print("h: ", h)
        if error < tol:
            w = w3
            t = t + h
            T.append(t)
            W.append(w)
            H.append(h)
            
            if error < tol/128:
                h = 2*h

            if t+h > b:
                h = b - t

        else:
            h = h/2
            if h < min_h:
                print("Minimum h exceeded, step size too small")
                return 
            
    return T, W, H


def Problem12B():
    a = 0
    b = 1600
    w = 8
    tol = 10E-4
    T, W, H = RK4Var(a, b, w, tol)

    N = len(T)
    #plot(T,W)
    

    print(f"Step Size has become too small at N = {len(T)}")
    
    
    return T, W, H


def Problem13A():
    a = 0
    b = 2
    w = 3
    tol = 10E-4
    T, W, H = RK4Var(a, b, w, tol)

    N = len(T)
    print(f"Step Size has become too small at N = {len(T)}")
     
    return T, W, H

if __name__=='__main__':
    
    # a = 0
    # b = 0.2
    # w = 0.1
    # min_h = 10E-4
    # tol = 10E-4

    # T,W,H = RK4Var(a, b, w, tol,
    #         h_init = 0.1,
    #         min_h = 10E-4, 
    #         some_function=chemical_function)

    # T, W, H = Problem12B()
    # T,W,H = Problem13A()

    # for t,w,h in zip(T,W,H):
    #     print(f"{t:.1f}\t{h:.1f}\t{w:.5f}")