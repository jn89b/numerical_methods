#include <iostream>
#include <euler_method.h>

//euler method main

float approximate_function(float x, float y)
{   
    //from the hw y'=y-x^2+1
    return y-(x*x)+1;
}

int main()
{
    const double h = 0.2;

    double y0 = 0.5;
    double x0 = 0.0;
    int N = 10; //number of iterations

    //test for loop for function c++ 
    for (int i = 0; i < N; i++)
    {
        std::cout << "y = " << y0 << std::endl;
        std::cout << "x = " << x0 << std::endl;
        std::cout << std::endl;

        y0 = y0 + h * approximate_function(x0, y0);
        x0 = x0 + h;

    }

    // // create an object of EulerMethod
    // EulerMethod euler(0, 1, 0.1, 0.5);
    // // compute the y
    // double y = euler.compute_y();
    // // print the result
    // std::cout << "y = " << y << std::endl;
    return 0;
}