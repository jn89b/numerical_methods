clear;
clc;
close all;


n = 9
a  = 1;
b = 3;
alpha = 17;
beta = (43/3);
h = (b-a)/(n+1)


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

tolerance = 10E-8;
while max_Fw > tolerance

    i = i+1;
    if i > 10
        disp('too many iterations')
        break
    end

    [c,d,e] = J(W,alpha,beta,h);
    P = tridiag(c,d,e,-Fw);
    W = W + P;
    max_Fw = max(compute_F(X,W,alpha,beta,h));
    MF = [MF; max_Fw];
end

function Fw = compute_F(X,W,alpha,beta,h)
    n = length(X);
    Fw = zeros(n,1);
    
    Fw(1) = -alpha + (2*W(1)) - W(2) + h^2*(4+1/4*X(1)^3-W(1)*W(2)-alpha/(16*h));
    for i = 2:n-1
        Fw(i) = -W(i-1) + (2*W(i)) - W(i+1) + h^2*(4+1/4*X(i)^3-W(i)*W(i+1));
    end
    Fw(n) = -W(n-1) + (2*W(n)) - beta + h^2*(4+1/4*X(n)^3-W(n)*beta+beta/(16*h));
    
endfunction

function [c,d,e] = J(W,alpha,beta,h)
    n = length(W);
    c = zeros(n,1);
    d = zeros(n,1);
    e = zeros(n,1);
    
    d(1) = 2 - (h/16) * (W(2)-alpha);
    e(1) = -1 - h * (W(1)/16);
    for i = 2:n-1
        c(i) = -1 + h * (W(i)/16);
        d(i) = 2 - (h/16) * (W(i+1)-W(i-1));
        e(i) = -1 - h * (W(i)/16);
    end

    c(n) = -1 + h * (W(n)/16);
    d(n) = 2 - (h/16) * (beta-W(n-1));
endfunction

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

endfunction

