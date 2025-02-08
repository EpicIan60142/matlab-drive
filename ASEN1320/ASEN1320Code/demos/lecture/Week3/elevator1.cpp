#include <iostream>

using namespace std;

int main()
{
  
  int    floor;

  // Ask user to enter a floor number
  cout << "Enter Floor : ";
  cin  >> floor;

  // Print out an error message if user enters 13, and exit 
  if (floor == 13 || floor == 666) 
  { 
    cout << "Error: There is no such floor like " << floor << endl;
    return 1;   //  returning with error
  }
  
  cout << "The elevator will travel to the floor " << floor << endl;

  return 0;     //  returning with no error
}