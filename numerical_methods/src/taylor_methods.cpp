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
    double best_h;
    double error_tolerance = 1e-4;
    double best_N;

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
    return compute_approx_dy(x, y) - 2*x - 2;
}

double compute_approx_ddddy(double x, double y)
{
    return compute_approx_dy(x, y) - 2*x - 2;
}

double compute_error(double actual_y, double y)
{
    return abs(y-actual_y);
}


double compute_taylor2(const double &h, double x, double y)
{
    return y + (h * compute_approx_dy(x, y)) 
        + ((pow(h,2) /2) * compute_approx_ddy(x, y));

}

double compute_taylor4(const double &h, double x, double y)
{
    return y + (h * compute_approx_dy(x, y)) 
        + ((pow(h,2) /2) * compute_approx_ddy(x, y))
        + ((pow(h,3) /6) * compute_approx_dddy(x, y))
        + ((pow(h,4) /24) * compute_approx_ddddy(x, y));
}

double compute_bound_error(double x, double factor_number=6.0)
{
    return 0.1 * (pow(EULER,2.0) - 2.0) * (pow(EULER, x) - 1.0) / factor_number;
}

double compute_bound_error_estimate(const double&h, double ME, double x, 
    double factor_number=6.0, double raise_number=2.0)
{
    return ME * (pow(EULER, x) - 1.0) * pow(h,raise_number) / (factor_number);
}
 
double compute_test(double x, double y)
{
    return y - pow(x,2) - (2*x) - 1;
}

double compute_true_d5y(double x_final=2.0)
{
    return 0.5*pow(EULER, x_final);
}

double compute_best_h(double error_tolerance, double ME, double factor_number=6.0, double order=4.0)
{
    return pow((error_tolerance)/((ME/factor_number) * (pow(EULER, 2.0) - 1.0)), order);
}


taylorInfo get_taylor2_history(const double &h, double &y0, double &x0, int N)
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
    double L = 1.0;

    for (int i=0; i<=N; i++)
    {
        taylorHistory.bounded_error(i) = abs(compute_bound_error(taylorHistory.x(i)));
        taylorHistory.estimated_bounded_error(i) = abs(compute_bound_error_estimate(h, 
            max_estimated_error, taylorHistory.x(i), 6));
    
    }

    return taylorHistory;
}

taylorInfo get_taylor4_history(const double &h, double &y0, double &x0, int N)
{
    taylorInfo taylorHistory = initialize_taylorinfo(h, y0, x0, N);
    Eigen::VectorXd estimated_error_margin = Eigen::VectorXd::Zero(N+1);

    for (int i =0; i<=N; i++)
    {
        taylorHistory.x(i) = x0;
        taylorHistory.y(i) = y0;
        taylorHistory.true_y(i) = actual_function(x0);
        taylorHistory.local_error(i) = compute_error(taylorHistory.true_y(i), y0);
        estimated_error_margin(i) = abs(compute_approx_ddddy(x0, y0));
        y0 = compute_taylor4(h, x0, y0);
        x0 += h;
    }

    //get max error margin 
    double max_estimated_error = estimated_error_margin.maxCoeff();
    double max_error = compute_true_d5y();

    for (int i=0; i<=N; i++)
    {
        taylorHistory.bounded_error(i) = abs(compute_bound_error_estimate(h, 
            max_error, taylorHistory.x(i), 120.0, 4.0));

        taylorHistory.estimated_bounded_error(i) = abs(compute_bound_error_estimate(h, 
            max_estimated_error, taylorHistory.x(i), 120.0, 4.0));

        taylorHistory.best_h = compute_best_h(taylorHistory.error_tolerance, 
            max_estimated_error, 120.0, 1/4.0);

        taylorHistory.best_N = ceil((taylorHistory.x(N) - taylorHistory.x(0)) / 
            taylorHistory.best_h);
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
        std::cout<< taylorHistory.x(i) << "\t" << 
        taylorHistory.y(i) << "\t" << taylorHistory.local_error(i) << "\t" <<
        taylorHistory.estimated_bounded_error(i) << std::endl;
    }
}

void taylor_4_5_3_A()
{
    const double h = 0.2;
    double y0 = 0.5;
    double x0 = 0.0;
    int N = 10; //number of iterations
    taylorInfo taylorInfo4 = get_taylor4_history(h, y0, x0, N);
    print_taylor(taylorInfo4, "taylor4");

}

void taylor_4_5_3_B(){
    const double h = 0.2;
    double y0 = 0.5;
    double x0 = 0.0;
    int N = 10; //number of iterations
    taylorInfo taylorInfo4 = get_taylor4_history(h, y0, x0, N);
    print_taylor(taylorInfo4, "taylor4");
}



Eigen::VectorXd get_other_element(const Eigen::VectorXd &vector, int element)
{
    int n = vector.size();
    int m = (n-1)/ element + 1;

    Eigen::VectorXd other_element(m);

    for (int i = 0; i < n; i+=element)
    {
        other_element(i/element) = vector(i);
    }

    return other_element;

}

void taylor_4_5_3_C(){
    const double h = 0.2;
    const double h2 = 0.1;
    const double h3 = 0.05;
    const double h4 = 0.025;
    
    double decimal_precision = 4;

    double y0 = 0.5;
    double x0 = 0.0;
    
    int N = 10; //number of iterations
    int N2 = 20;
    int N3 = 40;
    int N4 = 80;

    taylorInfo taylorHistory1 = get_taylor4_history(h, y0, x0, N);
    taylorInfo taylorHistory2 = get_taylor4_history(h2, y0, x0, N2);
    taylorInfo taylorHistory3 = get_taylor4_history(h3, y0, x0, N3);
    taylorInfo taylorHistory4 = get_taylor4_history(h4, y0, x0, N4);

    Eigen::VectorXd other_element2 = get_other_element(
        taylorHistory2.local_error, N2/10);

    Eigen::VectorXd other_element3 = get_other_element(
        taylorHistory3.local_error, N3/10);

    Eigen::VectorXd other_element4 = get_other_element(
        taylorHistory4.local_error, N4/10);


    std::cout << taylorHistory2.local_error << std::endl;

    //divide other_element2 by other_element3
    Eigen::VectorXd ratio1 = taylorHistory1.local_error.array() / other_element2.array(); 
    Eigen::VectorXd ratio2 = other_element2.array() / other_element3.array();
    Eigen::VectorXd ratio3 = other_element3.array() / other_element4.array();

    std::cout << std::setprecision(decimal_precision);
    std::cout<< "h" << "\t" << "R1" << "\t" << \ 
    "R2" << "\t"  << "R3" << std::endl;

    for (int i = 0; i <= N; i++)
    {
        //if i = 0 continue 
        if (i == 0) continue;

        std::cout<< taylorHistory1.x(i) << "\t"  << \
        ratio1(i) << "\t" << ratio2(i) << "\t" << ratio3(i) << std::endl;
    }

}

void taylor_4_5_3_D()
{
    //compute step size needed to get error less than error tolerance
    double error_tolerance = 1e-4;
    const double h = 0.2;
    double y0 = 0.5;
    double x0 = 0.0;
    int N = 10; //number of iterations

    taylorInfo taylorHistory = get_taylor4_history(h, y0, x0, N);    
    taylorInfo bestTaylor = get_taylor4_history(taylorHistory.best_h, 
        y0, x0, taylorHistory.best_N);

    //get maximum error in the best euler method
    double max_error = bestTaylor.local_error.maxCoeff();

    std::cout << "Number of steps N: " << bestTaylor.best_N << std::endl;
    std::cout << "optimal step size h: " << bestTaylor.best_h << std::endl;
    std::cout << "max error: " << max_error << std::endl;
   
}


int main()
{
    const double h = 0.1;
    double y0 = 0.5;
    double x0 = 0.0;
    int N = 20; //number of iterations

    taylorInfo taylorInfo4 = get_taylor4_history(h, y0, x0, N);
    print_taylor(taylorInfo4, "taylor4");print_taylor(taylorInfo4, "taylor4");
    taylor_4_5_3_A();
    std::cout << std::endl;

    taylor_4_5_3_B();
    std::cout << std::endl;

    taylor_4_5_3_C();
    std::cout << std::endl;
    
    taylor_4_5_3_D();

}