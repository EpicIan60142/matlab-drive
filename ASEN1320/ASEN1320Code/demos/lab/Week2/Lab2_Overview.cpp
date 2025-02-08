#include <iostream>
#include <string>

using namespace std;

void output();
void vars_and_types();
void operators();
void input();
void more_operations();

int main() {
	output();
	vars_and_types();
	operators();
	input();
	more_operations();
}

void output() {
    cout << "Welcome back to ASEN " << 1320 << "!" << "\n";
    cout << "This is Section 102" << endl; // can also use endl
    cout << "Let's start Recitation #2!\n\n";
}

void vars_and_types() {
	
	// declaring variables
	char c;		// c = ?
	int i;		// i = ?
	double d;	// d = ?
	string s;	// s = ?
	bool b;		// b = ?
	
	// undefined behavior if variables are used before declaration (can't assume initalized to zero)
	cout << "i = " << i << " and d = " << d << endl;
	
	// initializing variables
	c = 'a';							cout << "c = " << c << endl;
	i = 9;								cout << "i = " << i << endl;
	i = 5.5;							cout << "i = " << i << endl; // POLL: doubles casted to ints get truncated
	d = 1.2;							cout << "d = " << d << endl;
	s = "strings use DOUBLE quotes";	cout << "s = " << s << endl;
	b = true;							cout << "b = " << b << "\n\n"; // lowercase 't' in true and 'f' in false

	// special characters
	char newline = '\n';
	char tab = '\t'; 
	cout << "Here's a newline:" << newline << "and a tab:" << tab << "!\n\n";

	// constants
	const double pi = 3.14159;
	
	// ERRORS
	// double a; // redeclaring a variable
	// pi = 3; // redifing constants
}

void operators() {

	int a, b; // a = ? and b = ?
	
	a = 5;
	b = 2;
	cout << "a = " << a << " and b = " << b << endl;
	
	// operations
	cout << "a + b = " << a + b << endl;
	cout << "a - b = " << a - b << endl;
	cout << "a * b = " << a * b << endl;
	cout << "a / b = " << a / b << endl;					// POLL: dividing ints gets floored
	cout << "(double)a / b = " << (double)a / b << endl;	// casting
	cout << "a % b = " << a % b << "\n\n";

}

void input() {

	// input with strings
	string name;
    cout << "What's your name? ";
    cin >> name;
	cout << "Hello, " << name << "!\n\n";
	
	// input with numbers
	int birth_year;
	cout << "What is your birth year? ";
	cin >> birth_year;

	// operations using inputs
	int age;
	age = 2020 - birth_year;
	cout << "Wow! You're around " << age << " years old!\n";
	cout << "Wow! You're around " << 2020 - birth_year << " years old!\n\n"; // you can do operations inline
	
	// multiple inputs
	string s;
	int i;
	cout << "Input a word and a number: ";
	cin >> s >> i;
	cout << "s = \"" << s << "\" and i = " << i << "\n\n";
}

void more_operations() {
	int i = 1;
	bool result;
	
	// assignment operators
	cout << "i = " << i << endl;
	i++;		cout << "i = " << i << " after i++\n";
	i++;		cout << "i = " << i << " after i++\n";
	i--;		cout << "i = " << i << " after i--\n";
	i+=2;		cout << "i = " << i << " after i+=2\n";
	i-=3;		cout << "i = " << i << " after i-=3\n\n";
	
	// comparisons
	result = 1 == 1;	cout << "1 == 1 --> " << result << endl;
	result = 1 != 2;	cout << "1 != 2 --> " << result << endl;
	result = 1 > 1;		cout << "1 > 1  --> " << result << endl;
	result = 1 >= 1;	cout << "1 >= 1 --> " << result << endl;
	result = 1 < 2;		cout << "1 < 2  --> " << result << endl;
	result = 1 <= 2;	cout << "1 <= 2 --> " << result << "\n\n";
	
	// logical operators
	result = !true;				cout << "!true = " << result << endl;
	result = !false;			cout << "!false = " << result << endl;
	result = true && true;		cout << "true && true = " << result << endl;
	result = true && false;		cout << "true && false = " << result << endl;
	result = false && false;	cout << "false && false = " << result << endl;
	result = true || true;		cout << "true || true = " << result << endl;
	result = true || false;		cout << "true || false = " << result << endl;
	result = false || false;	cout << "false || false = " << result << endl;
}