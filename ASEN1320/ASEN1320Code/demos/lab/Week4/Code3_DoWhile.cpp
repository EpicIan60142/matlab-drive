#include <iostream>

using namespace std;

int main()
{
    // Print frist n natural number and its sum;
    
    int number;         // to store the number
    int sum;            // to store the sum
    
    cout<<"Enter a number: ";
    cin>>number;
    
    sum=0;
    
    // Do While loop from i=1 to i=number; 
    
    int i =1;
   do
    {
        cout<<i<<" ";       // Print i 
        sum += i;           // Add i to sum
        i++;                // increment i
    
    }
    while(i<=number);           // semicolon
    
    cout<<"\nSum of numbers from 1 to " << number <<" is " << sum<<endl;
    
    
    
    
    
}