#include "utilities.h"
// g++ -o run main.cpp utilities.cpp
// ./run
int main () {
   int const Number = 4;
   double angle = 90.0; // [deg] 
   double loc[Number] = {1, 0, 0, 1};
   rotate(loc, Number, angle);  
   write_csv(loc, Number, "output.csv"); 
   return 0;
}