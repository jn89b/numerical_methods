#include <iostream>

double function(double &x, double &y)
{   
    return y-(x*x)+1;
}

void euler(const double &h, double &y0, double &x0, int &N)
{
    //test for loop for function c++ 
    for (int i = 0; i < N; i++)
    {
        std::cout << "y = " << y0 << std::endl;
        std::cout << "x = " << x0 << std::endl;
        std::cout << std::endl;

        y0 = y0 + h * function(x0, y0);
        x0 = x0 + h;
    }
}

int main()
{
    const double h = 0.2;
    double y0 = 0.5;
    double x0 = 0.0;
    int N = 10; //number of iterations

    euler(h, y0, x0, N);

    return 0;
}