#include <fstream>
#include <iostream>
#include <cmath>
#include <string>
#include "utilities.h"

using namespace std;

// Write a function named "rotate" that takes 1D array of type double 
// of a certain length, and modify values as instructed   
void rotate(double loc[], int length, double theta){
    double PI = 3.141592653589793;
    int RowNumber = length/2;
    for (int i = 0; i < RowNumber; i++) {
      double xrot = cos(theta*PI/180.0) * loc[i*2] - sin(theta*PI/180.0) * loc[i*2+1];
      double yrot = sin(theta*PI/180.0) * loc[i*2] + cos(theta*PI/180.0) * loc[i*2+1];
      loc[i*2] = xrot; 
      loc[i*2+1] = yrot;   
    }
}

// Write a function named "write_csv" that writes out values stored in 1D array 
// of type double of a certain lenth to cvs file as instructed. 
void write_csv(double output[], int length, string filename){
    int RowNumber = length/2;
    ofstream outputFile;  
    outputFile.open(filename);
    for (int i = 0; i < RowNumber; i++) outputFile << output[i*2] << ", " << output[i*2+1] << endl;
    outputFile.close();
} 