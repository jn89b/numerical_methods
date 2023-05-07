clear;
clc;
close all;

n = 999;
a  = 1;
b = 3;
alpha = 17;
beta = (43/3);
h = (b-a)/(n+1);

X = zeros(n,1);
W = zeros(n,1);

for i = 1:n
    X(i) = a + i*h;
    W(i) = ((beta-alpha)/(b-a)) * (X(i)-a) + alpha;
end

i = 1;
Fw = compute_F(X,W,alpha,beta,h);
max_Fw = max(Fw);
MF = max_Fw;
tolerance = 10^-8;

while max_Fw > tolerance
    i = i+1;
    if i > 10
        disp('too many iterations')
        break
    end

    [c,d,e] = J(X,W,alpha,beta,h);
    P = tridiag(c,d,e,-Fw);
    W = W + P;
    Fw = compute_F(X,W,alpha,beta,h);
    max_Fw = max(abs(Fw));
    MF = [MF; max_Fw];
end
y = compute_true_function(X);
error = abs(y-W);
max_error_w = max(abs(y-W));
format shortg 
%print max error together 
disp('max error in W is')
disp(max_error_w)
x_t_w_y_error = [X W y error]

%% Functions used 
function true_val = compute_true_function(x)
    true_val = x.^2 + (16./x);
end

function f = compute_function(x,y,yprime)
    f = 4 + (1/4)*x^3 - y*yprime/8;
end

function fy = compute_dfunction(x,y,yprime)
    fy = -yprime/8;
end

function fyy = compute_ddfunction(x,y,yprime)
    fyy = -y/8;
end

function Fw = compute_F(X,W,alpha,beta,h)
    n = length(X);
    Fw = zeros(n,1);
    
    dy1 = (W(2)-alpha)/(2*h);
    Fw(1) = -alpha + (2*W(1)) - W(2) + h^2*compute_function(X(1),W(1), dy1);
    
    for i = 2:n-1
        dyi = (W(i+1)-W(i-1))/(2*h);
        Fw(i) = -W(i-1) + (2*W(i)) - W(i+1) + h^2*compute_function(X(i),W(i), dyi);
    end
    
    dyn = (beta-W(n-1))/(2*h);
    Fw(n) = -W(n-1) + (2*W(n)) - beta + h^2*compute_function(X(n),W(n), dyn);
    
end

function [c,d,e] = J(X,W,alpha,beta,h)
    n = length(W);
    c = zeros(n,1);
    d = zeros(n,1);
    e = zeros(n,1);
    
    d(1) = 2 + h^2 * compute_dfunction(X(1), W(1), (W(2)-alpha)/(2*h));
    e(1) = -1 + h^2 * compute_ddfunction(X(1), W(1), (W(2)-alpha)/(2*h)) * (1/(2*h));
    
    for i = 2:n-1
        c(i) = -1 + h^2 * compute_ddfunction(X(i), W(i), (W(i+1)-W(i-1))/(2*h)) * (-1/(2*h));
        d(i) = 2 + h^2 * compute_dfunction(X(i), W(i), (W(i+1)-W(i-1))/(2*h));
        e(i) = -1 + h^2 * compute_ddfunction(X(i), W(i), (W(i+1)-W(i-1))/(2*h)) * (1/(2*h));
    end

    c(n) = -1 + h^2 * compute_ddfunction(X(n), W(n), (beta-W(n-1))/(2*h)) * (-1/(2*h));
    d(n) = 2 + h^2 * compute_dfunction(X(n), W(n), (beta-W(n-1))/(2*h));
end

function W = tridiag(c,d,e,b)
    n = length(d);
    W = zeros(n,1);
    for k = 2:n
        mult = -c(k)/d(k-1);
        d(k) = mult*e(k-1) + d(k);
        b(k) = mult*b(k-1) + b(k);
    end
    W(n) = b(n)/d(n);
    for k = n-1:-1:1
        W(k) = (b(k)-e(k)*W(k+1))/d(k);
    end

end

