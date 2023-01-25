#ifndef __EULER_METHOD_H__
#define __EULER_METHOD_H__

// https://lefticus.gitbooks.io/cpp-best-practices/content/03-Style.html

class EulerMethod
{
    public:
        // create a constructor
        EulerMethod(double x0, double y0, double h, double x);

        // create a destructor
        ~EulerMethod();

        //set values for x0, y0, h, x
        void set_values(double x0, double y0, double h, double x);

        //method to compute the y
        double compute_y();

    private:
        // variables
        double x0;
        double y0;
        double h;
        double x;

};


#endif