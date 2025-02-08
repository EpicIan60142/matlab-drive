#include <iostream>

using namespace std;

int main()
{
  int    floor;
  // Ask user to enter a floor number
  cout << "Enter Floor : ";
  cin  >> floor;

switch (floor)
{
  case 13: 
    cout << "Error: There is no thirteenth floor" << endl;
    break;
  case 666:
    cout << "Error: There is no 666 floor" << endl;
    break;
  default:
    cout << "The elevator will travel to the floor " << floor << endl;
}    

  return 0;     //  returning with no error
}
