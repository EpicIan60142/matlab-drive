#include <iostream>

using namespace std;

// Function Declaration
// Different Types  - int,double, void, etc.

double areaCalc(int r);  // Tell function return type and argument type can be different variable types
// double areaCalc(double); // Tell them this is also a valid declaration (need to define the variable name in defintion)


const double pi=3.1416;   // Explain global variable and scope

int main()
{
    
    int r = 5;
    double area;
    
    area = areaCalc(r);
    cout<< "Area using func = "<< area;
    
    
    return 0;
}


double areaCalc(int r)
{
    double area = pi*r*r;
    
    return area;
}
