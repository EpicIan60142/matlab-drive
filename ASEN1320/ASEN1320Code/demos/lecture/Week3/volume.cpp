#include<iostream>

using namespace std;

int main()
{

  const double PI = 3.141592653589793; 
  double volume, radius;

  cout << "Enter a radius: ";
  cin  >> radius;
 
  volume =  4.0 / 3.0 * PI *radius*radius*radius;

  cout << "The volume is "  << int(volume+0.5) << "\n";
  
  return 0;
}