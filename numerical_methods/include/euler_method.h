#ifndef __EULER_METHOD_H__
#define __EULER_METHOD_H__

// https://lefticus.gitbooks.io/cpp-best-practices/content/03-Style.html

class EulerMethod
{
    public:
        // create a constructor
        EulerMethod(double x0, double y0, double h, int N);

        // create a destructor
        ~EulerMethod();

        //set values for x0, y0, h, N
        void set_values(double x0, double y0, double h, int N);

        //method to compute the y
        double invoke_function(double &x, double &y, 
            double (*approximate_function)(double, double));
        
        double approximate_function(double &x, double &y);

        double compute_bound_error(double x);

        //method actual function 
        double actual_function(double &x);

        void compute_y();

    private:
        // variables
        double x0;
        double y0;
        double h;
        int N;

};


#endif