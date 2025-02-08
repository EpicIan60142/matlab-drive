#include <cstdlib>
#include "point.h"
// >>g++ -o point main.cpp point.cpp 
// >>./point 

int main() {
 Point A[5]; // initialize an array of Point object
 srand (1);  // initialize the random number generator seed 
 for( int i=0; i<5; i++){
     double a = rand()*(1.0 / RAND_MAX); // generate random numbers between 0 and 1
     double b = rand()*(1.0 / RAND_MAX); 
	 A[i].setPoint(a, b); 
	 A[i].setFlag(a, b);
	 A[i].displayInfo();
 }
  return 0;
}