clear;
clc;
close all;

x = 0.0
y = 0.01

chem_value = chemical_function(x,y)


function chem_value = chemical_function(x,y)
  k = 6.22E-19;
  n_1 = 2E3
  n_2 = 2E3
  n_3 = 3E3 
  
  chem_value = k*(n_1-y/2)^2*(n_2-(y/2))^2*(n_3-(3*y/4))^3
  
endfunction