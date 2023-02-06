#include <iostream>     // std::cout
#include <cmath>        // std::abs

/*

HW 3
Example c in Chapter 5.2
Example d in chapter 5.2

Example a in Chapter 5.3 
Example b in Chapter 5.3
due next wednesday 
*/

const double EULER = 2.71828182845904523536;

double function(double &x, double &y)
{   
    return y-pow(x,2)+1;
}

double dfunction (double &x, double &y)
{
    return y - pow(x,2) -2*x + 1;
}

double ddfunction (double &x, double &y)
{
    return function(x, y) - 2*x - 2;
}

double dddfunction (double &x, double &y)
{
    return function(x, y) - 2*x - 2;
}

void taylor2(const double &h, double &y0, double &x0, int &N)
{
    //test for loop for function c++ 
    std::cout<<"taylor2" << std::endl;
    std::cout<< "x" << "\t" << "y guess" << std::endl;

    for (int i = 0; i <= N; i++)
    {

        //cout as table format 
        std::cout<< x0 << "\t" << y0 << std::endl;
        y0 = y0 + (h * function(x0, y0)) + ((pow(h,2) /2) * dfunction(x0, y0));
        x0 = x0 + h;
    }
}

void taylor4(const double &h, double &y0, double &x0, int &N)
{

    std::cout<<"taylor4" << std::endl;
    std::cout<< "x" << "\t" << "y guess" << std::endl;

    for (int i = 0; i <= N; i++)
    {

        //cout as table format 
        std::cout<< x0 << "\t" << y0 << std::endl;

        y0 = y0 + (h * function(x0, y0)) + \
            ((pow(h,2) /2) * dfunction(x0, y0)) + \
            ((pow(h,3) /6) * ddfunction(x0, y0)) + \
            ((pow(h,4) /24) * dddfunction(x0, y0));
        
        x0 = x0 + h;
    }

}

int main()
{

    const double h = 0.2;
    double y0 = 0.5;
    double x0 = 0.0;
    int N = 10; //number of iterations

    taylor2(h, y0, x0, N);
    
    std::cout << std::endl;
    y0 = 0.5;
    x0 = 0.0;
    taylor4(h, y0, x0, N);

}