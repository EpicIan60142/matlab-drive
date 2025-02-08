//warn about copying each other's code -- we can tell!
    
//tell them we went easy for the first assignment, but in the future they need to start early and follow directions exactly
        //show "guide to good grades" slide from Dr. Matsuo

#include <iostream>
#include <string> //emphasize that they need this one if using strings (not on this hw, but will use later)

using namespace std;

int main()
{
    //if statements = decision statements
        //Draw flow chart...
            // if(some condition)
                //task
        //mention how if the condition is false, everything in brackets gets skipped
    //in practice...
    int b = 10;
    
    if (b < 11)
    {
        cout << "The first condition" << endl;
    }
    
    if (b > 10)
    {
        cout << "The second condition" << endl;
    }
    
    //run again with b = 15 or somehting and talk through
        //--> POTENTIAL POLL: what will the output be if I use b = 15?
    
    //That's a little clunky...
        //Draw another flow chart!
        //     start
        //       |
        // (if condition)    --> yes -->    <task 1>
        //       |
        //       --> no --> <task 2>
        
    //practical example...
    int age = 10;
    if (age >= 21)
    {
        cout << "allowed" << endl;
    } else {
        cout << "not allowed" << endl;
    }
    
    //can make it even more complicated with "nested if statements" aka if else!
        //flow chart...
        //if(condition 1)
        //      task 1
        //else if(condition 2)
        //      task 2
        //else
        //      task 3
    //that else takes care of all the leftovers
    
    int temp = 50;
    if (temp > 80)
    {
        cout << "it's hot outside" << endl;
    } else if (temp > 70 && temp <= 80)
    {
        cout << "the weather is perfect" << endl;
    } else 
    {
        cout << "it's cold outside" << endl;
    }
    //note -- only ONE of those three options will run
    //note -- talk about && for "and" and || for "or"
    
    //Practical example:
    double dividend, divisor;
    cout << "enter dividend: ";
    cin >> dividend;
    cout << "enter divisor: ";
    cin >> divisor;
    
    if(divisor==0)
    {
        cout << "ERROR! There is a zero in the denominator";
        return 1;
    }
    
    double result = dividend/divisor;
    
    cout << "The result is " << result << endl;
    
    //run that with divisible number and with a 0 denom to demonstrate error (which is in their homework)
    
    //mention creating if statements based on booleans ( review the idea of booleans )
    
    //------------------------------------------------------------------------------------------
    //That is all they should need for their homework.  from here on is only if you have time...
    //------------------------------------------------------------------------------------------
    
    //SWITCH CASES
    int choice = -1;
    
    switch (choice) {
        case 1:
            cout<<"we are in case 1" << endl;
            break;
        case 2:
            cout<<"we are in case 2" << endl;
            break;
        default : 
            cout<<"invalid case" << endl;
    }
    //talk about the break, the default case
    //these are good for when you have a set number of known options
    
    
    
    
    
    //optional practical example
    char grade = 'A';
    
    switch (grade) {
      case 'A' :
         cout << "Excellent work!" << endl; 
         break;
      case 'B' :
         cout << "Good job" << endl;
         break;
      case 'C' :
         cout << "Average" << endl;
         break;
      case 'D' :
         cout << "Passing" << endl;
         break;
      case 'F' :
         cout << "Failed" << endl;
         break;
      default :
         cout << "Invalid grade" << endl;
   }
   cout << "Your grade is " << grade << endl;
    
    return 0;
}
