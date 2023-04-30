import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


true_function = lambda x : x**2 + (16/x)

def set_pandas_display_options() -> None:
    """Set pandas display options."""
    # Ref: https://stackoverflow.com/a/52432757/
    display = pd.options.display

    display.max_columns = 1000
    display.max_rows = 1000
    display.max_colwidth = 199
    display.width = 1000

def vector_function(t,w):
    f1 = lambda x,w : w[1]
    f2 = lambda x,w : 4 + (0.25*x**3) - ((1/8)*w[0]*w[1]) 
    f3 = lambda x,w : w[3]
    f4 = lambda x,w : -((1/8)*w[1]*w[2]) - ((1/8)*w[0]*w[3])   
    
    fv = np.zeros(len(w))
    fv[0] = f1(t,w)
    fv[1] = f2(t,w)
    fv[2] = f3(t,w)
    fv[3] = f4(t,w)
    return fv

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

    for i in range(N):
        s1 = function(T[i], W[i])
        s2 = function(T[i] + h_val/2, W[i] + h_val*s1/2)
        s3 = function(T[i] + h_val/2, W[i] + h_val*s2/2)
        s4 = function(T[i] + h_val, W[i] + h_val*s3)
        W[i+1] = W[i] + (h_val*(s1 + 2*s2 + 2*s3 + s4)/6)
        T[i+1] = T[i] + h_val
    
    return T,W

def AdaptiveRKVector(a,b,w,tol,condition=1):
    T = []
    W = []
    H = []
    t = a 
    w = w
    h = 0.1
    
    if a + h > b:
        h = b - a
    
    while t < b - 10E-8:
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

def Problem17a():
    """
    RK_VEC_VAR_11_1_A
    """
    a = 1 
    b = 3 
    alpha = 17
    beta = 43/35
    s_current  = (beta-alpha)/(b-a)
    tol = 10E-4
    w = np.array([alpha,s_current,0,1])
    T,W,H = AdaptiveRKVector(a,b,w,tol,condition=1)

    W = np.array(W)
    W1 = W[:,0]
    W3 = W[:,2]

    g = np.zeros(len(W1))
    s = np.zeros(len(W1))

    g_current = W1[-1] - beta
    g[0] = g_current
    s[0] = s_current

    i = 0 
    while abs(g_current) > 10E-8:
         
        if i > 10:
            print("The method has failed to converge")
            break
         
        s_current = s[i] - ((W1[-1] - beta)/ W3[-1]) 
        s[i+1] = s_current
        
        w = np.array([alpha,s_current,0,1])
        T,W,H = AdaptiveRKVector(a,b,w,tol,condition=1)
        W = np.array(W)
        W1 = W[:,0]
        W3 = W[:,2]

        g_current = W1[-1] - beta
        g[i+1] = g_current
        i = i + 1
        
    final_answer_dict = {
        's': s,
        'g': g,
    }

    true_list = []
    error_list = []
    
    for t,w in zip(T,W1):
        
        true_val = true_function(t)
        true_list.append(true_val)

        error = abs(true_val - w)
        error_list.append(error)

    final_answer_df = pd.DataFrame(final_answer_dict)
    print(final_answer_df)

    info_dict = {
        'T': T,
        'H': H,
        'W1': W1,
        'true': true_list,
        'error': error_list,
    }

    info_df = pd.DataFrame(info_dict)
    print(info_df)
    
    return final_answer_dict

def Problem17B():
    """
    RK_VEC_VAR_11_2_A
    """
    a = 1 
    b = 3 
    alpha = 17
    beta = 43/3
    s_current  = (beta-alpha)/(b-a)
    tol = 10E-5
    N = 20
    w = np.array([alpha,s_current,0,1])
    T,W = RK4Vector(a,b,w,N)

    W = np.array(W)
    W1 = W[:,0]
    W3 = W[:,2]

    g = np.zeros(len(W1))
    s = np.zeros(len(W1))

    g_current = W1[-1] - beta
    g[0] = g_current
    s[0] = s_current

    i = 0 
    while abs(g_current) > tol: 
        if i > 10:
            print("The method has failed to converge")
            break
         
        s_current = s[i] - (W1[-1] - beta)/ W3[-1] 
        s[i+1] = s_current
        
        w = np.array([alpha,s_current,0,1])
        T,W = RK4Vector(a,b,w,N)

        W = np.array(W)
        W1 = W[:,0]
        W3 = W[:,2]

        g_current = W1[-1] - beta
        g[i+1] = g_current
        i = i + 1
        
    answer_dict = {
        'g': g,
        's': s,
    }

    answer_df = pd.DataFrame(answer_dict)
    print("final answer",answer_df)

    true_list = []
    error_list = []
    
    for t,w in zip(T,W1):
        
        true_val = true_function(t)
        true_list.append(true_val)

        error = abs(true_val - w)
        error_list.append(error)

    info_dict = {
        'T': T,
        'W1': W1,
        'true': true_list,
        'error': error_list,
    }

    info_df = pd.DataFrame(info_dict)
    print(info_df)

    return answer_df

def Problem18A():
    """
    RK_VEC_VAR_11_2_A
    """
    a = 1 
    b = 3 
    alpha = 17
    beta = 43/3
    s_current  = (beta-alpha)/(b-a)
    tol = 10E-4
    w = np.array([alpha,s_current,0,1])
    T,W,H = AdaptiveRKVector(a,b,w,tol,condition=1)

    W = np.array(W)
    W1 = W[:,0]
    W3 = W[:,2]

    g = np.zeros(len(W1))
    s = np.zeros(len(W1))

    g_current = W1[-1] - beta
    g[0] = g_current
    s[0] = s_current

    i = 0 
    while abs(g_current) > 10E-8:
         
        if i > 10:
            print("The method has failed to converge")
            break
         
        s_current = -20 + (0.1) * (i-1)#s[i] - ((W1[-1] - beta)/ W3[-1]) 
        print(s_current)
        s[i+1] = s_current
        
        w = np.array([alpha,s_current,0,1])
        T,W,H = AdaptiveRKVector(a,b,w,tol,condition=1)
        W = np.array(W)
        W1 = W[:,0]
        W3 = W[:,2]

        g_current = W1[-1] - beta
        g[i+1] = g_current
        i = i + 1
        
    final_answer_dict = {
        's': s,
        'g': g,
    }

    true_list = []
    error_list = []
    
    for t,w in zip(T,W1):
        
        true_val = true_function(t)
        true_list.append(true_val)

        error = abs(true_val - w)
        error_list.append(error)

    final_answer_df = pd.DataFrame(final_answer_dict)
    print("final answer", final_answer_df)

    info_dict = {
        'T': T,
        'H': H,
        'W1': W1,
        'true': true_list,
        'error': error_list,
    }

    info_df = pd.DataFrame(info_dict)
    print(info_df)
    
    return final_answer_dict
    
if __name__=='__main__':
    
    set_pandas_display_options()
    seventeen_dict = Problem17a()
    print('\n')
    seventeen_dict_b = Problem17B()
    print('\n')
    # solution = Problem18A()
    #%% 
    """
    RK_VEC_VAR_11_2_A
    """
    a = 1 
    b = 3 
    alpha = 17
    beta = 43/3
    s_current  = -20#(beta-alpha)/(b-a)
    tol = 10E-4
    w = np.array([alpha,s_current,0,1])
    T,W,H = AdaptiveRKVector(a,b,w,tol,condition=1)


    W = np.array(W)
    W1 = W[:,0]
    W3 = W[:,2]

    n = 300
    g = np.zeros(n+1)
    s = np.zeros(n+1)
    # g = np.zeros(len(W1))
    # s = np.zeros(len(W1))

    g_current = W1[-1] - beta
    g[0] = g_current
    s[0] = s_current

    i = 0 
    while abs(g_current) > 10E-8:
         
        if i >= n:
            print("The method has failed to converge")
            break
         
        s_current = -20 + (0.1) * (i)#s[i] - ((W1[-1] - beta)/ W3[-1]) 
        # s_current = s[i] - ((W1[-1] - beta)/ W3[-1]) #-20 + (0.1) * (i-1)
        print(-20 + (0.1) * (i-1))#s[i] - ((W1[-1] - beta)/ W3[-1]) 

        s[i+1] = s_current 
        
        w = np.array([alpha,s_current,0,1])
        T,W,H = AdaptiveRKVector(a,b,w,tol,condition=1)
        W = np.array(W)
        W1 = W[:,0]
        W3 = W[:,2]

        g_current = W1[-1] - beta
        g[i+1] = g_current
        i = i + 1
        
    final_answer_dict = {
        's': s,
        'g': g,
    }

    true_list = []
    error_list = []
    
    for t,w in zip(T,W1):
        
        true_val = true_function(t)
        true_list.append(true_val)

        error = abs(true_val - w)
        error_list.append(error)

    final_answer_df = pd.DataFrame(final_answer_dict)
    print(final_answer_df)

    info_dict = {
        'T': T,
        'H': H,
        'W1': W1,
        'true': true_list,
        'error': error_list,
    }

    info_df = pd.DataFrame(info_dict)
    print(info_df)

    plt.close('all')
    plt.plot(s,g, '-')
    #plot horizontal line at y = 0
    plt.axhline(y=0, color='r', linestyle='-')
    plt.scatter(seventeen_dict_b['s'], seventeen_dict_b['g'])
    plt.xlabel('s')
    plt.show()
    
    s_list = []
    g_list = []
    for i in range(30):
        s_val = -20 + (0.1)*(i-1)
        g_val = (s_val)
        
        s_list.append(s_val)
        g_list.append(g_val)
        
        