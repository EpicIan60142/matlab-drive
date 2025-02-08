#include <iostream>

using namespace std;

int main()
{
    int n=1;
    
    
     // While is entry control do-while is exit control
     
    while(n<1)   // While checks the condition before execution of 1st iteration
    {
        cout<<"While: This text will not be printed";
    }
    
    do   // Do-while checks the condition after execution of 1st iteration
    {
        cout<<"Do-While: This text will be printed";
    }
    while(n<1);
    return 0;
    
}