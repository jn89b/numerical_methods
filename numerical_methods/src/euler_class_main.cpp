#include <iostream>
#include <euler_method.h>

int main()
{
    double h = 0.2;
    double y0 = 0.5;
    double x0 = 0.0;
    int N = 10; //number of iterations

    EulerMethod euler(x0, y0, h, N);
    
    euler.set_values(x0, y0, h, N);
    euler.compute_y();
}