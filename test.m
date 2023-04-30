 % Finite difference method 
 % Approximate the solution of y"=(-2/x)y'+(2/x^2)y+ sin(lnx)/x^2
 % for 1<=x<=2 with y(1)=1 and y(2)=2.
 
 p = @(x) (-2/x);  
 q = @(x) (2/x^2);
 r = @(x) (sin(log(x))/x^2);
 
c2 = -0.0392070132
c1 = 1.1392070132
true_function = @(x) (c1*x) + (c2/(x^2)) ... 
    - ((0.3)*sin(log(x))) ...
    - ((0.1)*cos(log(x))) 
 
 aa = 1; bb = 2; alpha = 1; beta = 2; n=9;      % h = (bb-aa)/(n+1);   h=0.1
  
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
 a(1) = 2+h^2*q(x);
 b(1) = -1+0.5*h*p(x);
 d(1) = -h^2*r(x)+(1+0.5*h*p(x))*alpha;
 m = n-1;
 
 for i = 2 : m 
   x = aa+i*h;
   a(i) = 2+h^2*q(x);
   b(i) = -1+0.5*h*p(x);
   c(i) = -1-0.5*h*p(x);
   d(i) = -h^2*r(x);
 end
 
 x = bb-h;
 a(n) = 2+h^2*q(x);
 c(n) = -1-0.5*h*p(x);
 d(n) = -h^2*r(x)+(1-0.5*h*p(x))*beta;
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
 
 