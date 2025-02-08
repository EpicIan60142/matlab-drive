#include <iostream>

using namespace std;

int main()
{
  int    floor = 1;

  while (floor <= 20)
  {
     if (floor == 13) cout << "There is no thirteenth floor" << endl;
     cout << "The elevator is passing the floor " << floor << endl;
     floor++;   
  }

  return 0;     //  returning with no error
}