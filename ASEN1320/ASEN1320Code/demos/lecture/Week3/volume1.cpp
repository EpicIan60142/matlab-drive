
#include<iostream>
#include<cmath>

using namespace std;

int main()
{

  double volume, radius;

  cout << "Enter a radius: ";
  cin  >> radius;

  volume =  4.0/3.0 * M_PI * pow(radius,3);
 
  cout << "The volume is "  << round(volume) << "\n";
 
  return 0;
}