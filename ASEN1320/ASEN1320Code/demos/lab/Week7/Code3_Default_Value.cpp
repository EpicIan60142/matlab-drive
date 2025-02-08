#include <iostream>

using namespace std;


// Default values
int addFn(int a = 10,int b = 5,int c = 6);
// int addFn(int a,int b = 5,int c = 6);            // valid declaration
// int addFn(int a,int b,int c = 6);                // valid declaration
// int add(int a,int b = 5,int c);                  // invalid declaration
// int add(int a = 2,int b = 5,int c);              // invalid declaration


int main()
{
    
    cout<<"Calling the function using 0 argument(s) = " << addFn()<<endl;
    cout<<"Calling the function using 1 argument(s) = " << addFn(2)<<endl;
    cout<<"Calling the function using 2 argument(s) = " << addFn(2,8)<<endl;
    cout<<"Calling the function using 3 argument(s) = " << addFn(2,8,5)<<endl;
    
}

//int addFn(int a = 10,int b = 5,int c = 6) // Invalid statement.
int addFn(int a ,int b,int c)
{
    return a+b+c;
}

// If there is time show an example of function overloading
// Area example

//double areaCalc(double r);
//double areaCalc(int l, int b);

