#include <iostream>     // std::cout
#include <cmath>        // std::abs
#include <Eigen/Dense>  // Eigen::VectorXd
#include <iomanip>

/*

HW 3
Example c in Chapter 5.2
Example d in chapter 5.2

Example a in Chapter 5.3 
Example b in Chapter 5.3
due next wednesday 
*/



const double EULER = 2.71828182845904523536;

typedef struct taylorInfo
{
    int N;
    Eigen::VectorXd x;
    Eigen::VectorXd y;
    Eigen::VectorXd local_error;
    Eigen::VectorXd true_y;
    Eigen::VectorXd bounded_error;
    Eigen::VectorXd estimated_bounded_error;
    const double h;
} taylorInfo;


taylorInfo initialize_taylorinfo(const double&h, double &y0, double &x0, int &N)
{
    //instiantiate eulerInfo with 0 vectors
    Eigen::VectorXd x = Eigen::VectorXd::Zero(N+1);
    Eigen::VectorXd y = Eigen::VectorXd::Zero(N+1);
    Eigen::VectorXd local_error = Eigen::VectorXd::Zero(N+1);
    Eigen::VectorXd true_y = Eigen::VectorXd::Zero(N+1);
    Eigen::VectorXd bounded_error = Eigen::VectorXd::Zero(N+1);
    Eigen::VectorXd estimated_bounded_error = Eigen::VectorXd::Zero(N+1);

    taylorInfo taylorHistory = {N, x, y, 
    local_error, true_y, 
    bounded_error, estimated_bounded_error, h};
    
    return taylorHistory;
}


double actual_function(double x)
{
    return (x+1)*(x+1)-0.5*exp(x);
}

double compute_approx_dy(double x, double y)
{   
    return y-(x*x)+1;
}

double compute_approx_ddy(double x, double y)
{
    return compute_approx_dy(x, y) - 2*x;
}

double compute_approx_dddy(double x, double y)
{
    return compute_approx_ddy(x, y) - 2*x - 2;
}

double compute_approx_ddddy(double x, double y)
{
    return compute_approx_ddy(x, y) - 2*x - 2;
}

double compute_error(double actual_y, double y)
{
    return abs(actual_y - y);
}


double compute_taylor2(const double &h, double x, double y)
{
    return y + (h * compute_approx_dy(x, y)) 
        + ((pow(h,2) /2) * compute_approx_ddy(x, y));

}

double compute_bound_error(double x)
{
    return 0.1 * (pow(EULER,2.0) - 2.0) * (pow(EULER, x) - 1.0) / 6.0;
}

double compute_bound_error_estimate(const double&h, double ME, double x)
{
    
    //get up to 2 decimal places of ME
    // ME = round(ME * 100) / 100;

    return h * ME * (pow(EULER, x) - 1.0) * pow(h,2) / (6 * 2);

    // double val_round = 1e4;
    // return round(bound_error * val_round)/ val_round;

} 

double compute_test(double x, double y){
    return y - pow(x,2) - (2*x) - 1;
}

taylorInfo get_taylor2_history(const double &h, double &y0, double &x0, int &N)
{
    taylorInfo taylorHistory = initialize_taylorinfo(h, y0, x0, N);
    Eigen::VectorXd estimated_error_margin = Eigen::VectorXd::Zero(N+1);

    for (int i = 0; i <= N; i++)
    {
        taylorHistory.x(i) = x0;
        taylorHistory.y(i) = y0;
        taylorHistory.true_y(i) = actual_function(x0);
        taylorHistory.local_error(i) = compute_error(taylorHistory.true_y(i), y0);
        estimated_error_margin(i) = abs(compute_test(x0, y0));
        y0 = compute_taylor2(h, x0, y0);
        x0 += h;
    }

    //get max error margin
    double max_estimated_error = estimated_error_margin.maxCoeff();
    // double max_error = compute_true_ddy();
    double L = 1.0;

    for (int i=0; i<=N; i++)
    {
        taylorHistory.estimated_bounded_error(i) = abs(compute_bound_error_estimate(h, 
            max_estimated_error, taylorHistory.x(i)));
    }

    return taylorHistory;
}


void print_taylor(taylorInfo &taylorHistory, std::string method)
{
    //print out the values
    double decimal_precision = 4.0;
    std::cout << std::setprecision(decimal_precision);

    std::cout<< method << std::endl;
    std::cout<< "x" << "\t" << "y guess" << "\t" << \
    "E" << "\t" << "EBE" <<  std::endl;
    
    for (int i = 0; i <= taylorHistory.N; i++)
    {
        std::cout<< taylorHistory.x(i) << "\t" << \
        taylorHistory.y(i) << "\t" << taylorHistory.local_error(i) << "\t" \ 
        << taylorHistory.estimated_bounded_error(i) << std::endl;
    }
}

int main()
{
    const double h = 0.2;
    double y0 = 0.5;
    double x0 = 0.0;
    int N = 10; //number of iterations

    taylorInfo taylorInfo2 = get_taylor2_history(h, y0, x0, N);
    print_taylor(taylorInfo2, "taylor2");
    // five_3_a(h, y0, x0, N);


}