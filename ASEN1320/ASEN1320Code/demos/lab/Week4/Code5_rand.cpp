#include <iostream>
#include <cstdlib>
#include <cmath>

using namespace std;


int main()
{
    //srand(2);       //  with and without srand. default srand(1)
    double randNum = rand();
    
    cout<<"Output from rand(): " << randNum << endl;
    
    // multiple rand() generate new numbers
    double randNum01 = rand() * (1.0 / RAND_MAX);
    cout<<"Between 0 and 1: " << randNum01 << endl;
    
    // Between 10 and 20
    int uppLim = 20;
    int lowLim = 10;
    
    double rn = rand() * (1.0 / RAND_MAX);
    int randNum10_20 = round(1 + (uppLim - lowLim) * rn);
  
    
    // another way
    int r = rand();
    int randNum10_20_2 = r % ((uppLim -lowLim) + 1) + lowLim;
    
    
    cout<<"Between 10 and 20: " << randNum10_20 << endl;
    cout<<"Between 10 and 20(2nd method): " << randNum10_20_2 << endl;
    
}

