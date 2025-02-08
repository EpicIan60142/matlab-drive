#include "utilities.h"

//g++ -c main.cpp utilities.cpp
//g++ -o run main.o utilities.o

int main () {
   int const Number = 2;
   double angle = 90.0; // [deg] 
   double loc[Number] = {0, 1};
   rotate(loc, Number, angle);  
   write_csv(loc, Number, "output.csv");
   return 0;
}
