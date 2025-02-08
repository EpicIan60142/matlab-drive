#include <iostream>

using namespace std;

int main() 
{
  int sum = 0;
  int n = 1;

    while (n <= 10) 
    {
      sum += n;  // sum = sum + n
      n++; 
      cout << sum << endl;
    }
    
  return 0;
}