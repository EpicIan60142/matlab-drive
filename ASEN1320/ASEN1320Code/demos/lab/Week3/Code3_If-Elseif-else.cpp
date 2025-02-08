#include <iostream>

using namespace std;


int main()
{
    int temp = 75;
    
    if (temp >= 80)  // if greater than 80
    {
        cout << "it's hot outside" << endl;
        return 0;
        
    } 
    
    else if (temp > 70 && temp <= 80)  // greater than 70 AND less than or equal to 80 And ->&& OR-> ||
    {
        cout << "the weather is perfect" << endl;
        return 0;
    }
    
    // else if
    //..
    //..
    
    else 
    {
        cout << "it's cold outside" << endl;
        return 0;
    }
    
    cout<<"Done running";
    
}





