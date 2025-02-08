#include <fstream>
#include <cmath>
#include <string>
#include "utilities.h"

using namespace std;

// Write a function that takes in 1D array and rotate it 
void rotate(double loc[], int length, double theta){
    double PI = 2*acos(0.0);
    double thetar = theta*PI/180.0;
    double xrot = cos(thetar) * loc[0] - sin(thetar) * loc[1];
    double yrot = sin(thetar) * loc[0] + cos(thetar) * loc[1];
    loc[0] = xrot; 
    loc[1] = yrot;   
}

// Write a function that write values 1D array of size 2 as cvs file
void write_csv(double output[], int length, string filename){
    ofstream outputFile;  
    outputFile.open(filename, ios::app);
    outputFile << output[0] << ", " << output[1] << endl;
    outputFile.close();
} 
