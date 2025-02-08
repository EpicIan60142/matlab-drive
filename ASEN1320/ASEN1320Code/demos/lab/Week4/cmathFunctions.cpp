#include <iostream>
#include <cmath>

using namespace std;

int main()
{
    
    // Power and Root
    
    int i;
    
    i = 8;          cout<< "Square of "<<i<<" is "<< pow(i,2) << endl;                      // pow(5,2) -> 5^2 = 25
    i = 8;          cout<< "Cube of "<<i<<" is "  << pow(i,3) << endl << endl;
    
    i = 4;          cout<< "Square Root of "<<i<<" is "  << sqrt(i) << endl;                // sqrt(4) -> √4 = 2
    i = 4;          cout<< "Square Root(using pow) of "<<i<<" is "  << pow(i,0.5) << endl;  // Also pow(4,0.5) -> √4 = 2
    i = 8;          cout<< "Cube Root of "<<i<<" is "  << cbrt(i) << endl << endl;          // cbrt(4) -> 3√8 = 2   
    
    
    // Rounding
    
    double j;
    
    j = 2.4;        cout<< "j rounded to integer " << round(j) << endl;  
    j = 2.5;        cout<< "j rounded to integer " << round(j) << endl << endl;
    
    // floor() always round the number down to nearest integer
    j = 2.4;        cout<< "floor(j) is " << floor(j) << endl;  
    j = 2.5;        cout<< "floor(j) is " << floor(j) << endl << endl;
    
    // floor() always round the number up to nearest integer
    j = 2.4;        cout<< "ceil(j) is " << ceil(j) << endl;  
    j = 2.5;        cout<< "ceil(j) is " << ceil(j) << endl << endl;
    
    
    // C++ Constants
    
    cout<< "C++ pi variable M_PI " << M_PI << endl;
    cout<< "C++ pi/2 variable M_PI " << M_PI_2 << endl << endl;
    
    // Trignometric functions
    // Always radians
    
    cout << "Sin of pi/2 is " << sin(M_PI_2) << endl;
    cout << "Cos of pi/2 is " << cos(M_PI_2) << endl;
    cout << "Tan of pi/2 is " << tan(M_PI_2) << endl;
    cout << "Tan of pi/4 is " << tan(M_PI_4) << endl << endl;
    
    cout << "aSin of 1 is " << asin(1) << endl;
    cout << "aCos of 0 is " << acos(0) << endl;
    cout << "aTan of 1 is " << atan(1) << endl;
}