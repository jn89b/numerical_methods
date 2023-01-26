#include <euler_method.h>
#include <iostream>


EulerMethod::EulerMethod(double x0, 
    double y0, double h, int N)
{
    this->x0 = x0;
    this->y0 = y0;
    this->h = h;
    this->N = N;
}

EulerMethod::~EulerMethod()
{
}

void EulerMethod::set_values(double x0, double y0, double h, int N)
{
    this->x0 = x0;
    this->y0 = y0;
    this->h = h;
    this->N = N;
}

double EulerMethod::approximate_function(double &x, double &y)
{   
    return y-(x*x)+1;
}

void EulerMethod::compute_y()
{
    for (int i =0; i < N; i++)
    {

        double y = y0 + h * approximate_function(x0, y0);
        double x = x0 + h;

        std::cout << "y = " << y0 << std::endl;
        std::cout << "x = " << x0 << std::endl;
        std::cout<< std::endl;

        set_values(x, y, h, N);       

    }

    return;
}
