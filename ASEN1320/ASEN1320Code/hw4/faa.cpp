#include <iostream>
#include <fstream>
#include <string>

using namespace std;

int main()
{
    //Declare 1 D array
    string airport[30];
    string airport_temp;
    int average[30], fy17[30], fy18[30];
    int i, j, temp;
    fstream inputStream;
    
    inputStream.open("AirTraffic.txt"); // connect the inputStream variable to a text file

    //reading in data line by line from "AirTraffic.txt" using loop    
    cout << "Average airport operation data before sorting\n";
    for(i = 0; i < 30; i++) 
    {  inputStream >> airport[i] >> average[i] >> fy17[i] >> fy18[i];
       cout << airport[i] << " " << average[i] << endl;
    }
    
    //sorting 1D arrays in ascending order accoding to average array values
    cout << "Average airport operation data after sorting\n";
    for(i = 0; i < 30; i++) 
    {
      for(j=i+1; j < 30; j++)
      {
         if (average[i] > average[j] )
         {
          temp       = average[i];
          average[i] = average[j];
          average[j] = temp;
          temp       = fy17[i];
          fy17[i]    = fy17[j];
          fy17[j]    = temp;
           temp      = fy18[i];
          fy18[i]    = fy18[j];
          fy18[j]    = temp;
          airport_temp = airport[i];
          airport[i]   = airport[j];
          airport[j]   = airport_temp; 
          
         }
      }
      cout << airport[i] << " " << average[i]  << endl;
    }
    
    // Select airports with both fy17 and fy18 higher than average
     for(i = 0; i < 30; i++) 
    { 
        if (fy17[i] > average[i] && fy18[i] > average[i]) 
        cout << airport[i] << " airport operation rose above average in FY17-18" << endl;
    }
    
    inputStream.close(); // Closing a text file
    
    return 0;   
}