#include <iostream>

using namespace std;

int main() 
{
  int sum = 0;
  int n;
  
    for (n = 1; n <= 10; n++) 
    {
      sum += n;  // sum = sum + n 
      cout << sum << endl;
    }
    
  return 0;
}