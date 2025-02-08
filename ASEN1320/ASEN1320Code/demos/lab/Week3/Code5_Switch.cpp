#include <iostream>

using namespace std;


int main()
{
    int choice = 2;

    
    switch (choice)
    {
        case 1 :    cout<<"We are in Case 1\n";
                    break;    
        
        case 2 :    cout<<"We are in Case 2\n";
                    break;
        
        default :   cout<<"Invalid Choice\n";
    }
    
    // Same implementation as If-Else
    if(choice==1)
    {
        cout<<"IF Case: 1\n";
    }
    else if (choice == 2)
    {
        
        cout<<"IF Case: 2\n";
        
    }
    //else if
    
    else
    {
        cout<<"IF Case: Invalid\n";
    }
    
    return 0;
}



