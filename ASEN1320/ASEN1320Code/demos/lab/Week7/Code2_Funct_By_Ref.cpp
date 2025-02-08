#include <iostream>

using namespace std;

// Go through the swap example first to show simple call by reference

// DemoCode/Week6/lecture

// Example of void return
// Tell them you can have mix of reference and call by value.
void areaCalc(int r, double &area);  


const double pi=3.1416;   // Explain global variable and scope

int main()
{
    
    int r = 5;
    double area;
    
    areaCalc(r,area);
    cout<< "Value of r = " <<r <<endl;  // to show r wont chnage
    cout<< "Area using func = "<< area <<endl;
    
    
    return 0;
}


// No return here
void areaCalc(int r, double &area)
{
    area = pi*r*r;
    //cout<<area<<endl;   // show the difference between area and &area
    //cout<<&area<<endl;
    //r = r+5;  show r wont change, since its not a reference;
}