#include <iostream>

using namespace std;


int main()
{
    int num1 = 1;
    int num2 = 2;
    
    if(num1 > 10)    // First if else
    {
        cout << "Num1 is greater than 10\n";
        
        if(num2 > 10)   // Second if else inside first if.
        {
            cout << "Num2 is greater than 10\n";
        }
        // else if
        else
        {
            cout<<"Num2 is less than or equal to 10\n";
        }
    }
    //else if
    else   // part of first if else
    {
        cout<<"Num1 is less than or equal to 10\n";
        
        if(num2>10)   // Second if else inside first else
        {
            cout << "Num2 is greater than 10\n";
        }
        else
        {
            cout<<"Num2 is less than or equal to 10\n";
        }
    }
    
    
}


