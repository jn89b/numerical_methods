import numpy as np

"""
Problem RK4C 
"""

some_function = lambda x, y: y - x**2 + 1
actual_function = lambda x: (x+1)*(x+1) - 0.5*m.exp(x)


k = 6.22E-19
n_1 = 2E3
n_2 = 2E3
n_3 = 3E3 

chemical_function = lambda x,y:\
    k*(n_1-y/2)**2*(n_2-(y/2))**2*(n_3-(3*y/4))**3
    
r = 0.1 #ft
g = 32.1 #ft/s^2
water_flow_function = lambda x,y:\
    (-0.6*np.pi*r**2)*(np.sqrt(2*g)) * ((np.sqrt(y))/ (np.pi*y**2))  

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

def Problem9A():
    a = 0
    b = 0.2
    y0 = 0
    x0 = 0
    N = 10

    x, y, actual_vals = RK4(a, b, N, x0, y0, chemical_function, chemical_function)
    #round to nearest integer
    y_estimate_init = round(y[-1])
    
    #make adjustments to N
    N = N * 2
    x, y, actual_vals = RK4(a, b, N, x0, y0, chemical_function, chemical_function)
    y_adjusted_estimated = round(y[-1])

    print("y_estimate = ", y_estimate_init)
    while y_adjusted_estimated != y_estimate_init:
        y_estimate_init = y_adjusted_estimated
        print("y_estimate = ", y_estimate_init)
        N = N * 2
        x, y, actual_vals = RK4(a, b, N, x0, y0, chemical_function, chemical_function)
        y_adjusted_estimated = round(y[N]) 

    print("Best N = ", N)

def Problem9B():
    a = 0 
    b = 1600 #seconds 
    y0 = 8 
    N = 80

    x, y, actual_vals = RK4(a, b, N, a, y0, 
        water_flow_function, water_flow_function, set_complex=True)

    #get last 6 indices 
    x = x[-6:]
    y = y[-6:]
    actual_vals = actual_vals[-6:]

    print("h\t x\t y\t ")
    for i in range(len(x)):
        #print("x = ", x[i], "y = ", y[i], "actual = ", actual_vals[i])
        print("{:.2f}\t {:.2f}\t {:.2f}".format((b-a)/N, x[i], y[i]))

if __name__ == "__main__":

    #Problem9A()
    Problem9B()

