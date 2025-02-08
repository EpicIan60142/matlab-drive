#include <iostream> 
#include <fstream>
#include <string>

using namespace std;

int main() {

 ifstream File;
 File.open("file.csv");
 string stringX, stringY;
 
 for( int i=0; i<10; i++){
	 getline(File, stringX, ',' );
     getline(File, stringY );
	 double a = stod(stringX);
	 double b = stod(stringY);
	 cout << "X: " << a  << " Y: " << b << endl;
 }
	
  return 0;
}