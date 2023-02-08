#include <iostream>
//import eigen library
#include <Eigen/Dense>

//create struct for euler method
typedef struct eulerInfo
{
    int N;
    Eigen::VectorXd x;
    Eigen::VectorXd y;
    Eigen::VectorXd local_error;
    Eigen::VectorXd true_y;
    Eigen::VectorXd bounded_error;
    const double h;

} eulerInfo;

double actual_function(double &x)
{
    return (x+1)*(x+1)-0.5*exp(x);
}


double function(double &x, double &y)
{   
    return y-(x*x)+1;
}

double dfunction (double &x, double &y)
{
    return y - pow(x,2) -2*x + 1;
}


double compute_error(double &x, double &y)
{
    return std::abs(actual_function(x) - y);
}

eulerInfo euler(const double &h, double &y0, double &x0, int &N)
{
    Eigen::VectorXd y(N);
    Eigen::VectorXd x(N);
    Eigen::VectorXd local_error(N);
    Eigen::VectorXd true_y(N);
    Eigen::VectorXd bounded_error(N);

    double y_guess = y0;
    double x_current = x0; 

    //test for loop for function c++ 
    for (int i = 0; i < N; i++)
    {
        y_guess = y_guess + h * function(x_current, y_guess);
        x_current = x_current + h;
        x(i) = x_current;
        y(i) = y_guess;
        local_error(i) = compute_error(x_current, y_guess);
        true_y(i) = actual_function(x_current);
        bounded_error(i) = abs(dfunction(x_current, y_guess));
    }

    eulerInfo eulerHistory = {N, x, y, 
        local_error, true_y, bounded_error, h};

    return eulerHistory;
}

void print_euler(eulerInfo &eulerHistory, std::string method)
{
    //print out the values
    std::cout<< method << std::endl;
    std::cout<< "x" << "\t" << "y guess" << "\t" << \
    "local error" << std::endl;
    
    for (int i = 0; i < eulerHistory.N; i++)
    {
        std::cout<< eulerHistory.x(i) << "\t" << \
        eulerHistory.y(i) << "\t" << eulerHistory.local_error(i) \ 
        << std::endl;
    }
}

Eigen::VectorXd get_other_element(const Eigen::VectorXd &vector, int element)
{
    int n = vector.size();
    int m = (n-1)/ element + 1;

    Eigen::VectorXd other_element(m);

    for (int i = 0; i < n; i+=element-1)
    {

        other_element(i/element) = vector(i);
    }

    return other_element;

}



void five_2_c()
{
    const double h = 0.2;
    const double h2 = h/2;
    const double h3 = h2/2;
    const double h4 = h3/2;

    double y0 = 0.5;
    double x0 = 0.0;
    
    int N = 10; //number of iterations
    int N2 = 20;
    int N3 = 40;
    int N4 = 80;

    eulerInfo eulerHistory1 = euler(h, y0, x0, N);
    eulerInfo eulerHistory2 = euler(h2, y0, x0, N2);
    eulerInfo eulerHistory3 = euler(h3, y0, x0, N3);
    eulerInfo eulerHistory4 = euler(h4, y0, x0, N4);

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


    std::cout<< "h" << "\t" << "R1" << "\t" << "\t" << \ 
    "R2" << "\t" <<  "\t" << "R3" << std::endl;

    for (int i = 0; i < 10; i++)
    {
        std::cout<< eulerHistory1.x(i) << "\t"  << \
        ratio1(i) << "\t" << ratio2(i) << "\t" << ratio3(i) << std::endl;
    }

}

void five_2_d(){

    double error = 100;
    double error_tolerance = 1e-4;

    double h = 0.2;
    double y0 = 0.5;
    double x0 = 0.0;
    int N = 10; //number of iterations

    std::cout<< "h" << "\t"  << "\t" << "N" \ 
    "\t" << "max error" << "\t" << std::endl;

    eulerInfo eulerHistory1 = euler(h, y0, x0, N);
    while (error > error_tolerance)
    {

        eulerInfo eulerHistory = euler(h, y0, x0, N);
        error = eulerHistory.local_error.maxCoeff();

        // double ME; 
        double ME = eulerHistory.bounded_error.maxCoeff(&ME);

        h = (error_tolerance * 2)/ (ME * exp(2) - 1);
        
        N = ceil((2.0)/h);

        std::cout<< h << "\t" << N << "\t" << error << "\t" << std::endl;

    }   
}

int main()
{
    const double h = 0.2;
    double y0 = 0.5;
    double x0 = 0.0;
    int N = 10; //number of iterations


    five_2_c();
    // five_2_d();


    return 0;
}