#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include <iomanip>

const double EULER = 2.71828182845904523536;

typedef struct eulerInfo
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

} eulerInfo;


eulerInfo initialize_euler(const double &h, double &y0, double &x0, int &N)
{
    //instiantiate eulerInfo with 0 vectors
    Eigen::VectorXd x = Eigen::VectorXd::Zero(N+1);
    Eigen::VectorXd y = Eigen::VectorXd::Zero(N+1);
    Eigen::VectorXd local_error = Eigen::VectorXd::Zero(N+1);
    Eigen::VectorXd true_y = Eigen::VectorXd::Zero(N+1);
    Eigen::VectorXd bounded_error = Eigen::VectorXd::Zero(N+1);
    Eigen::VectorXd estimated_bounded_error = Eigen::VectorXd::Zero(N+1);

    eulerInfo eulerHistory = {N, x, y, 
    local_error, true_y, 
    bounded_error, estimated_bounded_error, h};
    
    return eulerHistory;

}

double compute_true_ddy()
{
    return 0.5*pow(EULER, 2)-2;
}

double compute_approx_ddy(double x, double y)
{
    return y - pow(x,2) -2*x + 1;
}

double compute_dy(double x, double y)
{   
    return y-(x*x)+1;
}

double compute_euler(const double&h, double x0, double y0)
{
    return y0 + h*compute_dy(x0, y0);
}

double actual_function(double x)
{
    return (x+1)*(x+1)-0.5*exp(x);
}

double compute_error(double y_true, double y_guess)
{
    return abs(y_true - y_guess);
}


double compute_bound_error(double x)
{
    return 0.1 * (0.5 * pow(EULER,2.0) - 2.0) * (pow(EULER, x) - 1.0);
}

double compute_bound_error_estimate(const double&h, double L, double ME, double x)
{
    return 0.1 * ME * (pow(EULER, x) - 1.0);
}

double compute_best_h(double error_tolerance, double ME)
{
    return (error_tolerance)/((ME/2) * (pow(EULER, 2.0) - 1.0));
}

eulerInfo get_euler_method_info(const double &h, double y0, double x0, int N)
{
    eulerInfo eulerHistory = initialize_euler(h, y0, x0, N);
    
    Eigen::VectorXd estimated_error_margin = Eigen::VectorXd::Zero(N+1);

    for (int i = 0; i <= N; i++)
    {
        eulerHistory.x(i) = x0;
        eulerHistory.y(i) = y0;
        eulerHistory.true_y(i) = actual_function(eulerHistory.x(i));
        eulerHistory.local_error(i) = compute_error(eulerHistory.true_y(i), eulerHistory.y(i));
        estimated_error_margin(i) = abs(compute_approx_ddy(x0, y0));
        y0 = compute_euler(h, eulerHistory.x(i), eulerHistory.y(i)); 
        x0 = eulerHistory.x(i) + h;
    }

    //get max error margin
    double max_estimated_error = estimated_error_margin.maxCoeff();
    double max_error = compute_true_ddy();
    double L = 1.0;

    for (int i = 0; i <= N; i++)
    {
        eulerHistory.bounded_error(i) = compute_bound_error(eulerHistory.x(i));
        eulerHistory.estimated_bounded_error(i) = abs(compute_bound_error_estimate(
            h, L, max_estimated_error, eulerHistory.x(i))
        );
        
        eulerHistory.best_h = compute_best_h(eulerHistory.error_tolerance, 
            max_estimated_error);
        eulerHistory.best_N = ceil((eulerHistory.x(N) - eulerHistory.x(0))/
            eulerHistory.best_h);
    }

    return eulerHistory;
} 

void print_euler_info(eulerInfo eulerHistory, const int decimal_precision)
{
    std::cout << std::setprecision(decimal_precision);
    
    int N = eulerHistory.N;

    std::cout << "h" << "\t" << "x" << "\t" << "y" << "\t" 
    << "true y" << "\t" << "E" << "\t" << "EB" 
    << "\t" << "EBE" << std::endl;

    for (int i=0; i<=N; i++){
        std::cout << eulerHistory.h << "\t" << eulerHistory.x(i) << 
        "\t" << eulerHistory.y(i) << "\t" 
        << eulerHistory.true_y(i) << "\t" << eulerHistory.local_error(i) << "\t" << 
        eulerHistory.bounded_error(i) << "\t" << 
        eulerHistory.estimated_bounded_error(i) << std::endl;
    }
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

void five_2_B(){
    const double h = 0.2;
    double y0 = 0.5;
    double x0 = 0.0;
    int N = 10; //number of iterations

    eulerInfo eulerHistory = get_euler_method_info(h, y0, x0, N);
    print_euler_info(eulerHistory, 4);
}

void five_2_C(){
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

    eulerInfo eulerHistory1 = get_euler_method_info(h, y0, x0, N);
    eulerInfo eulerHistory2 = get_euler_method_info(h2, y0, x0, N2);
    eulerInfo eulerHistory3 = get_euler_method_info(h3, y0, x0, N3);
    eulerInfo eulerHistory4 = get_euler_method_info(h4, y0, x0, N4);

    Eigen::VectorXd other_element2 = get_other_element(
        eulerHistory2.local_error, N2/10);

    Eigen::VectorXd other_element3 = get_other_element(
        eulerHistory3.local_error, N3/10);

    Eigen::VectorXd other_element4 = get_other_element(
        eulerHistory4.local_error, N4/10);

    //divide other_element2 by other_element3
    Eigen::VectorXd ratio1 = other_element2.array() / eulerHistory1.local_error.array();
    Eigen::VectorXd ratio2 = other_element3.array() / other_element2.array();
    Eigen::VectorXd ratio3 = other_element4.array() / other_element3.array();

    std::cout << std::setprecision(decimal_precision);
    std::cout<< "h" << "\t" << "R1" << "\t" << \ 
    "R2" << "\t"  << "R3" << std::endl;

    for (int i = 0; i <= N; i++)
    {
        //if i = 0 continue 
        if (i == 0) continue;

        std::cout<< eulerHistory1.x(i) << "\t"  << \
        ratio1(i) << "\t" << ratio2(i) << "\t" << ratio3(i) << std::endl;
    }

}

void five_2_D(){
    //compute step size needed to get error less than error tolerance
    double error_tolerance = 1e-4;
    const double h = 0.2;
    double y0 = 0.5;
    double x0 = 0.0;
    int N = 10; //number of iterations

    eulerInfo eulerHistory = get_euler_method_info(h, y0, x0, N);    
    eulerInfo bestEuler = get_euler_method_info(eulerHistory.best_h, y0, x0, 
        eulerHistory.best_N);

    //get maximum error in the best euler method
    double max_error = bestEuler.local_error.maxCoeff();

    std::cout << "Number of steps N: " << eulerHistory.best_N << std::endl;
    std::cout << "optimal step size h: " << eulerHistory.best_h << std::endl;
    std::cout << "max error: " << max_error << std::endl;

}

int main(){
    
    five_2_B();
    std::cout << std::endl;
    five_2_C();
    std::cout << std::endl;
    five_2_D();


    return 0;
}


