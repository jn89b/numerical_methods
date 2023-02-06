#include <euler_method.h>
#include <iostream>     // std::cout
#include <cmath>        // std::abs

const double EULER = 2.71828182845904523536;

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

double EulerMethod::actual_function(double &x)
{
    return (x+1)*(x+1)-0.5*exp(x);
}

double EulerMethod::invoke_function(double &x, double &y, 
    double (*approximate_function)(double, double))
{
    return approximate_function(x, y);
}

double EulerMethod::compute_bound_error(double x)
{
    return 0.1 * (0.5 * pow(EULER,2.0) - 2.0) * (pow(EULER, x) - 1.0);
}

void EulerMethod::compute_y()
{
    for (int i =0; i < N; i++)
    {

        double y = y0 + h * approximate_function(x0, y0);
        double x = x0 + h;

        double actual_y = actual_function(x0);
        double diff =std::abs(actual_y - y0);
        double bounded_error = compute_bound_error(x0);

        std::cout << "w_i = " << y0 << std::endl;
        std::cout << "x = " << x0 << std::endl;
        std::cout << "actual y = " << actual_y << std::endl;
        std::cout << "error = " << diff << std::endl;
        std::cout << "bounded error = " << bounded_error << std::endl;
        std::cout<< std::endl;

        set_values(x, y, h, N);       

    }

    return;
}
