#include <iostream>

using namespace std;

int main()
{
// Variables declaration and Initialization//	   
   double I, V, R, C = 320.00;
   double eff;

// Ask for data entry for I and V from console
   cout << "Enter a resistance value: " ;
   cin  >> R;
   cout << "Enter a current value: " ;
   cin  >> I;
   cout << "Enter a voltage value: " ;
   cin  >> V;

// User input verification. If failed, return 
   if (I == 0.0 || V == 0.0 || R == 0.0)  
  {
    cout << "Error: resistance, current and voltage values cannot be zero!" << endl;
    return 1; // return with error
  }

/* Branching according to the maximum efficiency condition
   If the maximum efficiency condition is met, compute the maximum efficiency 
   Output respective messages to cosole (std output) */
   if (I*I > C/R)
   {
      cout << "The load is too high\n" ;
   }
   else if (I*I == C/R)
   {
      cout << "The maximum efficiency condition is met\n";
      eff =  (V*I - 2*I*I*R) / (V*I) * 100.00;                    // compute the efficiency 
      cout << "The efficiency is " << int(eff + 0.50) << "%\n" ; // 0.50 is added for rounding before type casting
   }
   else
   {
      cout << "The load is too low\n";
   }

   return 0; // return with no error
}