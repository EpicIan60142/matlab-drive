#include <iostream>     // include to use cout and cin
#include <fstream>      // include to use file streams
#include <string>       // include to use string

using namespace std;

int main() {

    /*********************************/
    /*            ARRAYS             */
    /*********************************/
    
    // What is an array?
    // --> a series of elements of the same type with a location relative to each other
    // --> it's still a variable so you must declare it with a name and a type
    // --> only now you must specify a length
    
    // declaring an array with name my_array of type int and of length 5
    int my_array[5] = {2, 4, 6, 8, 10};
    
    // illegal initialization!
    //my_array = {6, 7, 8, 9, 10};
    
    // ------------------------
    // ACCESSING values --> draw out the boxes and index numbers on the iPad!!!
    cout << "my_array: " << my_array[0] << " " << my_array[1] << " " << my_array[2] << " " << my_array[3] << " " << my_array[4] << endl;
    
    // in a loop!
    int i;
    cout << "my_array: ";
    for (i = 0; i < 5; i++) {
        cout << my_array[i] << " ";
    }
    cout << endl;
     
    // ------------------------
    // CHANGING VALUES
    my_array[4] = 7;     // POLL - what will the new list look like? {2, 7, 6, 8, 10} or {2, 4, 6, 7, 10} or {2, 4, 6, 8, 7}
    
    // let's check! --> again draw out the boxes
    cout << "my_array: ";
    for (i = 0; i < 5; i++) {
        cout << my_array[i] << " ";
    }
    cout << endl << endl;    
    
    // ------------------------
    // declared but not initalized --> sometimes you don't know what the values are initially!
    int empty_array[2];
    
    // remember that you CANNOT use a value before it is initalized
    cout << "empty_array: ";
    for (i = 0; i < 2; i++) {
        cout << empty_array[i] << " ";
    }
    cout << endl;
    
    // ------------------------
    // initializing the long way
    empty_array[0] = 8;
    empty_array[1] = 18;
    
    // now we see the expected values
    cout << "empty_array: ";
    for (i = 0; i < 2; i++) {
        cout << empty_array[i] << " ";
    }
    cout << endl << endl;
    
    
    // ------------------------
    // Some fun array operations
    // foo[0] = a;
    // foo[a] = 75;
    // b = foo [a+2];
    // foo[foo[a]] = foo[2] + 5;
    
    /*********************************/
    /*            FILE I/O           */
    /*********************************/
    
    // how to access a file?
    // --> through a "stream" object
    // --> you've already used streams! (cout and cin)
    // --> cout is an OUTPUT stream that you WRITE TO - you send data to it
    // --> cin is an INPUT strem that you READ FROM - you get data from it
    // --> files can be both INPUT and OUTPUT streams
    
    // CREATE a file stream object (variable)
    ofstream output_only_stream;            // type to output to a file
    ifstream input_only_stream;             // type to input from a file
    fstream input_and_output_stream;        // type to output to and input from a file
    
    // let's make a well-named stream for file input
    ifstream input_file;
    
    // ------------------------
    // OPEN a file stream
    input_file.open("input_file.txt");
    
    // check that the file is open!
    if (!input_file.is_open()) {
      cout << "Failed to open input_file.txt" << endl;
      return 1;
    }
    
    // ------------------------
    // READ from the file - same arrows as READING from cin!
    string word1, word2, word3;
    int num;
    
    // one word at a time with " " (space) as a delimiter (separator)
    input_file >> word1;
    input_file >> word2;
    input_file >> word3;
    input_file >> num;
    cout << word1 << " " << word2 << " " << word3 << " " << num << endl;
    
    // all at the same time
    input_file >> word1 >> word2 >>  word3 >> num;
    cout << word1 << " " << word2 << " " << word3 << " " << num << endl;
    // POLL: what is this going to print?
    //       "Welcome to ASEN 1320" or "Today is September 22"
    // answer: it will continue reading from where it was before!
    
    // ------------------------
    // WRITE to a file - same arrows as WRITING to cout!
    ofstream output_file;
    output_file.open("output_file.txt");
    if (!output_file.is_open()) {
      cout << "Failed to open output_file.txt" << endl;
      return 1;
    }
    
    // NOTE THE ENDL!!!!
    output_file << "Adding new text to this file!" << endl;
    
    // ------------------------
    // CLOSE a file stream - IMPORTANT!
    input_file.close();
    output_file.close();


    /*********************************/
    /*     EVERYTHING TOGETHER       */
    /*********************************/
    // I'm planning to build on this same code to demo sorting next week probably
    // The "Grades.txt" is like a gradebook for a student with a category, overll score, and percent of the overall grade

    // arrays to store data
    int n_grades = 5;
    string category[n_grades];
    double score[n_grades];
    int percent[n_grades];
    
    // input stream for reading from fi;e
    ifstream grades_file;
    
    // open file for reading
    grades_file.open("Grades.txt");
    if (!grades_file.is_open()) {
      cout << "Failed to open Grades.txt" << endl;
      return 1;
    }
    
    // read and store data line by line
    cout << endl << "The file contents:" << endl;
    for (int i = 0; i < n_grades; i++) {
      grades_file >> category[i] >> score[i] >> percent[i];
      cout << category[i] << " " << score[i] << " " << percent[i] << endl;
    }
    
    // calculate the student's final grade
    double final_grade = 0; // make sure to initialize your sum before using it
    for (int i = 0; i < n_grades; i++) {
      final_grade += score[i] * percent[i]/100.0;
    }
    cout << endl << "The final grade is " << final_grade << endl;
    
    // close file
    grades_file.close();
    
    return 0;
}
