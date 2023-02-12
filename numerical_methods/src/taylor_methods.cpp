#include <iostream>     // std::cout
#include <cmath>        // std::abs
#include <Eigen/Dense>  // Eigen::VectorXd
/*

HW 3
Example c in Chapter 5.2
Example d in chapter 5.2

Example a in Chapter 5.3 
Example b in Chapter 5.3
due next wednesday 
*/

typedef struct taylorInfo
{
    int N;
    const double h;
    double x0;
    double y0;
    std::string method;
    Eigen::VectorXd x;
    Eigen::VectorXd y;
    Eigen::VectorXd local_error;
    Eigen::VectorXd true_y;
    Eigen::VectorXd bounded_error;
} taylorInfo;


const double EULER = 2.71828182845904523536;

double actual_function(double &x)
{
    return (x+1)*(x+1)-0.5*exp(x);
}

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

double compute_error(double &x, double &y)
{
    return std::abs(actual_function(x) - y);
}

// double compute_bound_error(double &x, const double &h)
// {
//     return (h/2) * (0.5 * pow(EULER,2.0) - 2.0) * (pow(EULER, x) - 1.0);
// }

double compute_bound_error(double M2, double &x, const double &h)
{
    return (M2/2) * (pow(EULER, x) - 1.0) * pow(h,2);
}

taylorInfo taylor2(const double &h, double &y0, double &x0, int &N)
{

    //store taylor info 
    Eigen::VectorXd y(N+1);
    Eigen::VectorXd x(N+1);
    Eigen::VectorXd local_error(N+1);
    Eigen::VectorXd true_y(N+1);
    Eigen::VectorXd bounded_error(N+1);

    Eigen::VectorXd ME_vector(N+1);
    
    double error = 0.0;
    double y_guess = y0;
    double x_current = x0;


    for (int i = 0; i <= N; i++)
    {

        x(i) = x_current;
        y(i) = y_guess;
        local_error(i) = error;
        true_y(i) = actual_function(x_current);
        // bounded_error(i) = compute_bound_error(x_current, h);

        error = compute_error(x_current, y_guess);
        y_guess = y_guess + (h * function(x_current, y_guess)) + ((pow(h,2) /2) * dfunction(x_current, y_guess));
        x_current = x_current + h;

        ME_vector(i) = ddfunction(x_current, y_guess);

    }

    //get last value of x vector
    double M2 = ddfunction(x(N), y(N));
    

    for (int i = 0; i <= N; i++)
    {
        bounded_error(i) = abs(compute_bound_error(M2, x(i), h));      
        
    }


    taylorInfo taylor2History = {N, h, x0, y0, "taylor2", x, y, 
        local_error, true_y, bounded_error};

    std::cout<< "returning taylor2" << std::endl;

    return taylor2History;

}


void taylor4(const double &h, double &y0, double &x0, int &N)
{

    double error = 0.0;
    std::cout<<"taylor4" << std::endl;
    std::cout<< "x" << "\t" << "y guess" << "\t" << \
    "error" << std::endl;
    
    for (int i = 0; i <= N; i++)
    {
        //cout as table format 
        std::cout<< x0 << "\t" << y0 << 
        "\t" << error << std::endl; 
        
        error = compute_error(x0, y0);

        y0 = y0 + (h * function(x0, y0)) + \
            ((pow(h,2) /2) * dfunction(x0, y0)) + \
            ((pow(h,3) /6) * ddfunction(x0, y0)) + \
            ((pow(h,4) /24) * dddfunction(x0, y0));
        
        x0 = x0 + h;   
    }
}


void print_taylor(taylorInfo &taylorHistory, std::string method)
{
    //print out the values
    std::cout<< method << std::endl;
    std::cout<< "x" << "\t" << "y guess" << "\t" << \
    "local error" << "\t" << "bounded error" <<  std::endl;
    
    for (int i = 0; i <= taylorHistory.N; i++)
    {
        std::cout<< taylorHistory.x(i) << "\t" << \
        taylorHistory.y(i) << "\t" << taylorHistory.local_error(i) << "\t" \ 
        << taylorHistory.bounded_error(i) << std::endl;
    }
}

void five_3_a(const double &h, double &y0, double &x0, int &N)
{
    double error = 0.0;
    std::cout<<"five_3_b" << std::endl;
    
    // std::cout<< "x" << "\t" << "y guess" << "\t" << \
    // "error" << std::endl;
    taylorInfo taylorInfo2 = taylor2(h, y0, x0, N);
    print_taylor(taylorInfo2, "taylor2");

}

int main()
{
    const double h = 0.2;
    double y0 = 0.5;
    double x0 = 0.0;
    int N = 10; //number of iterations

    five_3_a(h, y0, x0, N);

    // taylorInfo taylorInfo2 = taylor2(h, y0, x0, N);

    // print_taylor(taylorInfo2, "Hello");

}