#include <euler_method.h>

EulerMethod::EulerMethod(double x0, 
    double y0, double h, double x)
{
    this->x0 = x0;
    this->y0 = y0;
    this->h = h;
    this->x = x;
}

EulerMethod::~EulerMethod()
{
}

void EulerMethod::set_values(double x0, 
    double y0, double h, double x)
{
    this->x0 = x0;
    this->y0 = y0;
    this->h = h;
    this->x = x;
}


double EulerMethod::compute_y()
{
    double y = this->y0;
    double x = this->x0;
    while (x < this->x)
    {
        y = y + this->h * (x + y);
        x = x + this->h;
    }
    return y;
}
