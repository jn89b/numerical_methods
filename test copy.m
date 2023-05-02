 % Finite difference method 
 % Approximate the solution of y"=(-2/x)y'+(2/x^2)y+ sin(lnx)/x^2
 % for 1<=x<=2 with y(1)=1 and y(2)=2.
 

 function [w] = f(x,y,yprime)
  w = 4 + 0.25*x^3 - y*yprime/8
end


function [w] = fy(x,y,yprime)
  w = -yprime/8
end


function [w] = fyprime(x,y,yprime)
  w = -y/8
end

true_function = @(x) x^2 + (16/x);

aa = 1; 
bb = 3; 
alpha = 17; 
beta = 43/35; 
n=29;      
  
 fprintf('   x           w   \n');
 h = (bb-aa)/(n+1);
 a = zeros(1,n+1);
 b = zeros(1,n+1);
 c = zeros(1,n+1);
 d = zeros(1,n+1);
 l = zeros(1,n+1);
 u = zeros(1,n+1);
 z = zeros(1,n+1);
 w = zeros(1,n+1);
 x = aa+h;

 w(1) = aa;
 w(n+1) = bb;

 for i = 1 : n 
  w(i) = aa + (i*h*(bb-aa))/(bb - x);
  x = aa+i*h;
 end

t = (w(1)-aa)/(2*h);
 a(1) = 2+h^2*fy(x,w(1),t);
 b(1) = -1+0.5*h*fyprime(x,w(1),t);
 d(1) = -((2*w(1))-w(2)-aa+(h^2*f(x,w(1),t)));
 m = n-1;
 
 for i = 2 : m 
   x = aa+i*h;
   a(i) = 2+h^2*fy(x,w(i),t);
   b(i) = -1+0.5*h*fyprime(x,w(1),t);
   c(i) = -1-0.5*h*fyprime(x,w(1),t);
   d(i) = -((2*w(i))-w(i-1)-w(i+1)+(h^2*f(x,w(i),t)))
 end
 
 x = bb-h;
 a(n) = 2 + h^2*fy(x,w(n),t);
 c(n) = -1-0.5*h*fyprime(x,w(n),t);
 d(n) = -((2*w(n))-w(n-1)-bb+(h^2*f(x,w(n),t)));
 l(1) = a(1);
 u(1) = b(1)/a(1);
 z(1) = d(1)/l(1);
 
 for i = 2 : m 
   l(i) = a(i)-c(i)*u(i-1);
   u(i) = b(i)/l(i);
   z(i) = (d(i)-c(i)*z(i-1))/l(i);
 end
 
 l(n) = a(n)-c(n)*u(n-1);
 z(n) = (d(n)-c(n)*z(n-1))/l(n);
 w(n) = z(n);

 for j = 1 : m 
   i = n-j;
   w(i) = z(i)-u(i)*w(i+1);
 end
 i = 0;
 fprintf('%5.4f    %11.8f\n', aa, alpha);
 for i = 1 : n
   x = aa+i*h;
   true_val = true_function(x);
   error = abs(true_val - w(i));
   fprintf('%5.4f    %11.8f %11.8f  11.8f ' , x, w(i), true_val, error);
   true
 end
 i = n+1;
 fprintf('%5.4f    %11.8f\n', bb, beta);
 



 