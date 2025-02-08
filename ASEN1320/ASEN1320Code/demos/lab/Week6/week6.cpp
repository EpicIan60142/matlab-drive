//good job on interview grading everybody!
//homework due on FRIDAY, same time as quiz

#include <iostream>
#include <fstream>
#include <string>

using namespace std;

//void change(int &n); //put this in later, only if you end up with time for a pass by reference example

int main()
{
    //Some people were confused about reading data from a text file from last week
    //Grading Example
    ifstream gradesinput; //initialize
    gradesinput.open("grades.txt"); //open
    if(!gradesinput.is_open()){ //check to make sure it opened successsfully
        cout<<"failed to open file!"<<endl;
        return 1;
    }
    
    //we're going to set values into arrays
    //initialize the arrays
        //look at grades file, there will be 5 of each, empty for now
    string sectionName[5];
    double myGrade[5];
    int classWeight[5];
    
    //read data in line-by-line with a loop
    cout << "read in the data:" << endl;
    for(int i=0;i<5;i++){
        gradesinput >> sectionName[i] >> myGrade[i] >> classWeight[i];
        cout << sectionName[i] << " " << myGrade[i] << " " << classWeight[i] << endl;
    }
    
    //use this to display one specific value if they seem confused about the arrays in the loop:
        //or for a potential poll question!
    //cout << "\ntest value: " << sectionName[3] << endl; 
    
    //we can do calculations with the arrays too! (now that we've read in data and set array values equal to the data)
    int totalExamWeights = classWeight[3] + classWeight[4];
    cout << "midterm and final exams account for " << totalExamWeights << " percent of final grade" << endl;
    
    
    //in the homework, you will read in data line by line just like that
    //and then you will sort the data using something called BUBBLE SORT
        //gone over in class, here is another example
    
    //very simple Bubble Sort example:
    int x[8] = { 9, 1, 4, 2, 5, 8, 4, 7 }; //initialize
    int j, k; // initialize increments
    int temp;
    
    cout << "Before sorting: ";
    for(int i=0; i<7; i++) cout<<" "<<x[i];
    cout << endl;
    
    for(j=0; j<7; j++){
        for(int k=j+1; k<7; k++){ //make sure you don't initialize in here! it will reset values every time - which we don't want
            // draw out the nested for loops! table with j,k,x[j],x[k] 
            if( x[j] > x[k] ){
                temp = x[j];
                x[j] = x[k];
                x[k] = temp;
            }
        }
    }
    
   cout << "After sorting: ";
   for(int i=0; i<7; i++) cout<<" "<<x[i];
   cout << endl;
   
   //in homework, you will read in data from a file, then run it through a bubble sort, combining the two things we just went over
        //can do this with the grades file we read in above!
        //go ahead and take the time to go over this in more detail if they seem confused
    
    
    
    
    
    
        
    //the rest of this isn't needed for homework, only do if you have time and students seem to understand how to bubble sort and/or read in data
    
    //call by reference
        //when we call by VALUE - makes a copy of waht is being passed
        //calling by reference passes actual memory address
        //so the function will have the acual n value, not just a number
    //simple example
        //go up to top and initialize: void change(int &n);
    //int n,sum;
    //cout << "enter the number n: " << endl;
    //cin >> n;
    //
    //change(n); //when function is called,n will become a reference to the argument
    //
    //Since a reference to a variable is treated exactly the same as the variable itself, changes made to the reference are passes through to the argument
    //cout << "Value of n: "<< n << endl;
    
    return 0;
}

//void change(int &n){
//    n=25;
//}










//On recitation topics: it would be helpful if you could 
//(1) have them practice use read in from a text file; 
//(2) set up a fixed size array and assign values in the loop; 
//(3) explain the fundamental difference between declaration and assignment - maybe via polls.  
    //They cannot put declaration in the loops for this HW and get OK on HW.
