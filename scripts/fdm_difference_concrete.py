import numpy as np
import pandas as pd

a = 1 
b = 3

alpha = 17
beta = 43/3
n = 9


x = np.zeros(n)
h = (b-a)/(n+1)
for i in range(n):
    x[i] = a + i * (b-a)/(n+1)

W = np.zeros(n)
#find linear line between a and b
for i in range(n):
    W[i] = alpha + (beta - alpha)/(b-a) * (x[i] - a)