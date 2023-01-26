#include <iostream>
#include <euler_method.h>


double approximate_function(double &x, double &y)
{   
    //from the hw y'=y-x^2+1
    return y-(x*x)+1;
}


int main()
{
    double h = 0.2;
    double y0 = 0.5;
    double x0 = 0.0;
    int N = 10; //number of iterations

    EulerMethod euler(x0, y0, h, N);
    euler.set_values(x0, y0, h, N);
    euler.compute_y();
    // euler.compute_y( &approximate_function(x0,y0));

}