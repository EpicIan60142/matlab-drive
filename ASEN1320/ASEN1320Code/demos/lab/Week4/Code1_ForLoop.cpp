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
    
    // For loop from i=1 to i=number; 
    // Show i can be decalred outside
    // Show you can have i+=2/3/4 also
    
    for(int i=1;i<=number;i++)
    {
        cout<<i<<" ";       // Print i 
        sum += i;           // Add i to sum
    
    }
    
    cout<<"\nSum of numbers from 1 to " << number <<" is " << sum<<endl;
    
    
    
    
    
}