#include "point.h"

void Point::setPoint(double x, double y) {
	xpoint = x;
	ypoint = y;
}
//Assign True/False if/not x is larger than y 
void Point::setFlag(double x, double y) {
  flagpoint = x > y; 
}
//Display private member data values
void Point::displayInfo() {
  std::cout << "x: "    << xpoint << "\n";
  std::cout << "y: "    << ypoint << "\n";
  std::cout << "flag: " << flagpoint << "\n";
}