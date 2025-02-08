#include <iostream>
#include <string> 

using namespace std;

int main()
{	
// Variables Definition & Constatnt variable definition and assignment//	
	const double PI = 3.141592653589793; 	
	string username;
	double density, radius, volume;
	int density_record, radius_record, volume_record;
	
// Greetings from HAL as cosole IO//
	cout << "Hello, my name is HAL!\nWhat is your name?\n";
	cin  >> username; 
	cout << "Hello, " << username << ". I am glad to meet you!\n";
	
// Data entry from standard input//	
	cout << "Enter a radius: ";
	cin  >> radius;
	cout << "Enter a density: ";
	cin  >> density;
	
// Volume computation and assignment //
    volume = 4.0 / 3.0 * PI *radius*radius*radius;
	
/* Type conversion and assignment  
   Note that 0.50 is added for rounding up and down to the nearest integer*/	
    density_record = density + 0.50;           // implicit type casting
    radius_record  = int(radius + 0.50);       // (old fashioned) explicit type casting 
    volume_record  = static_cast<int>(volume); // explicit type casting
  
// Output stream to standard ouput//    
    cout << "The radius is "  << radius_record  << "\n";
    cout << "The density is " << density_record << "\n";
    cout << "The volume is "  << volume_record  << "\n";
    
	return 0; 
}