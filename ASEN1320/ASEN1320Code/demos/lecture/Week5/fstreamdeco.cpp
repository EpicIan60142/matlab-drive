#include <iostream>
#include <fstream>

using namespace std;

int main()
{
	int number;
	ifstream inputExample;
	inputExample.open("Number.txt");
	inputExample >> number;
	cout << "ASEN " << number << " Hello, World!\n";
	inputExample.close();
	return 0; 
}