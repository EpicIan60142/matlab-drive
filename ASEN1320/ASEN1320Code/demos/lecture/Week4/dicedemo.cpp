#include<iostream>
#include<cstdlib>
#include<cmath>

using namespace std;

int main ()
{        
    int i, r;
    double rn; 
    for (i = 1; i <= 10; i++)
    {
       rn = rand() * (1.0 / RAND_MAX);  // 0-1
       r  = round(1 + (6 - 1) * rn); 
       cout <<  r << endl;
    }
    return 0; 
}