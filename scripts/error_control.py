import numpy as np
# import RKHomework89 
import math as m
# RK4 = RKHomework89.RK4

"""
These are assignments HW 11a and 11b 

"""

"""
HW 11A
Outputs are 
H           E1              E2               R
0.1         1.6595E-7       1.0594E-8        15.666
0.05        5.1974E-9       3.2826E-10       15.833
0.025       1.6259E-10      1.0215E-11       15.917
"""

"""
HW 11B 
Outputs are 
H           OE1              OE2            OE3               R
0.1         0.01659       0.0010594     0.00023609         15.666
0.05        0.016632      0.0010504   0.00023377        15.833
0.025       0.016649      0.001046    0.00023238        15.917
"""

"""
HW 12A - 27 in 5.4
Molecule equations but use adaptive step size
"""

some_function = lambda x, y: y - x**2 + 1
actual_function = lambda x,y: (x+1)*(x+1) - 0.5*m.exp(x)


k = 6.22E-19
n_1 = 2E3
n_2 = 2E3
n_3 = 3E3 
chemical_function = lambda x,y:\
    k*(n_1-y/2)**2*(n_2-(y/2))**2*(n_3-(3*y/4))**3


def RK4(a:float, b:float, N:int, x0:float, y0:float, 
    actual_function=some_function,
    function=some_function,
    set_complex=False, 
    set_h=(False,0.1)):
    
    """Computes RK4 for a given function"""
    if set_complex == True:
        y = np.zeros(N+1, dtype=complex)
        actual_vals = np.zeros(N+1, dtype=complex)
    else:
        y = np.zeros(N+1)
        actual_vals = np.zeros(N+1)
    

    x = np.zeros(N+1)

    x[0] = x0
    y[0] = y0
    actual_vals[0] = actual_function(x0,y0)
    
    #round to ceiling of h
    if set_h[0] == False: 
        h = (b-a)/N
    else:
        h = set_h[1]

    for i in range(N):
        x[i+1] = x[i] + h
        s1 = function(x[i], y[i])
        s2 = function(x[i] + h/2, y[i] + (h*s1)/2)
        s3 = function(x[i] + h/2, y[i] + (h*s2)/2)
        s4 = function(x[i] + h, y[i] + h*s3)
        y[i+1] = y[i] + (h*(s1 + 2*s2 + 2*s3 + s4)/6)
        
        actual_vals[i+1] = actual_function(x[i+1], y[i+1])
    return x, y, actual_vals


def RK_Test1():
    H = np.zeros(3)
    E1 = np.zeros(3)
    E2 = np.zeros(3)
    R = np.zeros(3)

    h = 0.1 #step size
    set_h = (True, h)
    #N = 3
    a = 0 
    b = 0.5
    x0 = 0
    y0 = 0.5
    for i in range(3):
        set_h = (True, h)
        
        x1,w1,actual_value1 = RK4(a, b, 1, x0, y0,
                            function=some_function, 
                            actual_function=actual_function,
                            set_h=set_h) #RK4(0, 0.5, 1, h)
        
        x2,w2,actual_value2 = RK4(a, b, 2, x0, y0,
                            function=some_function,
                            actual_function=actual_function, 
                            set_h=(True, h/2)) #RK4(0, 0.5, 2, h/2) #take a half step for the solution
        
        H[i] = h #store the step size
        # print("h = ", h)
        yv = actual_value1[-1] #actual value store based on step size
        E1[i] = abs(yv - w1[-1]) #error
        E2[i] = abs(yv - w2[-1]) #error
        #R[i] = E1[i]/E2[i] #ratio
        R[i] = (yv- w1[-1])/(yv - w2[-1]) #ratio
        h = h/2 #decrease step size

    return H, E1, E2, R

def RK_Test2():
    """
    (w3 - w1) / (w3 - w2) = 16

    so we solve for w3
    w3 = (16*w2 - w1)/15 

    error wise:
    w1 is C1(h^5)
    w2 is C2(h^5)
    w3 is C3(h^6)    

    Divide by h to get same C values 
    """
    H = np.zeros(3)
    OE1 = np.zeros(3)
    OE2 = np.zeros(3)
    OE3 = np.zeros(3)
    R = np.zeros(3)

    h = 0.1 #step size
    set_h = (True, h)
    #N = 3
    a = 0 
    b = 0.5
    x0 = 0
    y0 = 0.5
    w3 = []

    for i in range(3):

        x1,w1,actual_value1 = RK4(a, b, 1, x0, y0,
                            function=some_function, 
                            actual_function=actual_function,
                            set_h=(True, h)) #RK4(0, 0.5, 1, h)
        
        x2,w2,actual_value2 = RK4(a, b, 2, x0, y0,
                            function=some_function,
                            actual_function=actual_function, 
                            set_h=(True, h/2)) #RK4(0, 0.5, 2, h/2) #take a half step for the solution
        
        print("w1 = ", w1)
        w3 = (16*w2[-1] - w1[-1])/15 #RK4(0, 0.5, 3, h/4) #take a quarter step for the solution
        H[i] = h #store the step size
        yv = actual_value1[-1] #actual value store based on step size
        OE1[i] = abs(yv - w1[-1])/h**5#error
        OE2[i] = abs(yv - w2[-1])/h**5 #error
        OE3[i] = abs(yv - w3)/h**6 #error

        #R[i] = E1[i]/E2[i] #ratio
        R[i] = (yv- w1[-1])/(yv - w2[-1]) #ratio
        h = h/2 #decrease step size

    return H, OE1, OE2,OE3, R

def AdaptiveRK(a,b,w, tol,
               h_init = 0.1,
               min_h = 10E-4,
               some_function = chemical_function):
    T = []
    W = []
    H = []
    t = a
    h = h_init
    min_h = min_h

    if a + h > b:
        h = b - a

    bound_tolerance = b - (10E-12)
    while t < (bound_tolerance):

        x1,w1,actual_value1 = RK4(a, b, 1, t, w,
                            function=chemical_function, 
                            actual_function=chemical_function,
                            set_h=(True, h)) #RK4(0, 0.5, 1, h)
 

        x2,w2,actual_value2 = RK4(a, b, 2, t, w,
                            function=chemical_function,
                            actual_function=chemical_function, 
                            set_h=(True, h/2))
        
        w3 = (16*w2[-1]- w1[-1])/15
        error = abs(w3 - w1[-1])

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
                return T, W, H
            
    return T, W, H

if __name__ == "__main__":
    H, E1, E2, R = RK_Test1()
    
    H, OE1, OE2, OE3, R = RK_Test2()

    print("H \t", "E1 \t", "E2 \t", "R \t")


    for i in range(3):
        decimal_place = 15
        
        h = H[i].round(decimal_place)
        e1 = E1[i].round(decimal_place)
        e2 = E2[i].round(decimal_place)
        r = R[i].round(decimal_place)
        print(h,"\t", e1, "\t", e2, "\t", r)


    # for i in range(3):
    #     print("H = ", H[i], "E1 = ", E1[i], "E2 = ", E2[i], "R = ", R[i])

    # print("\n\n")
    
    # H, OE1, OE2, OE3, R = RK_Test2()
    # for i in range(3):
    #     print("H = ", H[i], "OE1 = ", OE1[i], "OE2 = ", OE2[i], "OE3 = ", OE3[i], "R = ", R[i])


    T, W, H = AdaptiveRK(a=0.0, 
                         b=0.2, 
                         w=0, 
                         tol=0.1, 
                         h_init = 0.1, 
                         min_h = 10E-12)
    
    print("T \t", "W \t", "H")

    for t,w,h in zip(T,W,H):
        
        #round decimal places t,w,h
        decimal_places = 8
        t = round(t, decimal_place)
        w = round(w, decimal_place)
        h = round(h, decimal_place)

        print(t, "\t", w, "\t", h)

#plot
import matplotlib.pyplot as plt

plt.close('all')
plt.plot(T,W)
plt.show()