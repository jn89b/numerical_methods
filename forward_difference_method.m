clear;
clc;
close all;

format shortg 

al = 1; 
h = 0.1;
v = 0.0005;
t = 0.5;

FDM_12_2_A(al,h, v, t)
FDM_12_2_A(1,0.1,0.01,0.5)

function FDM_12_2_A(al, h, v, t)
    m = 1/h;
    iterations = t/v;
    s = al^2*v/h^2;
    x = zeros(m-1,1);
    
    %go horizontal
    for i =1:m-1
        x(i) = i*h;
    end
    
    w = sin(pi*x);
    wn = zeros(m-1,1);
    
    %go vertical
    A = compute_A(m-1, al, v, h);
    
    for j = 1:iterations
        wn = A * w;
        w = wn;
    end
    
    %compute true solution
    u = true_solution(x,t);
    e = abs(u-w);
    x_w_u_e = [x,w,u,e]    

end

function [u] = true_solution(x,t)
    u = exp(-pi^2*t)*sin(pi*x);
end

function [A] = compute_A(iters, alpha, t, h)
    A = zeros(iters);
    
    lambda = alpha^2*(t/h^2);

    A(1,1) = (1 - 2*lambda);
    A(1,2) = lambda;

    for k = 2:iters-1
        A(k,k+1) = lambda;
        A(k,k) = (1 - 2*lambda);
        A(k,k-1) = lambda;
        continue;
    end

    A(iters,iters-1) = lambda;
    A(iters,iters) = (1 - 2*lambda);

end
