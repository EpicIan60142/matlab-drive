//Welcome back to Cloud9!

#include <iostream>
#include <string>
#include <fstream>

using namespace std;

//write this part when instructed...(skip the first time through)
class ThanksgivingFood{
    public: //can be accessout outside class
        string Name;
        int quantity;
        void displayInfo(){
            cout << "Food Name: " << Name << endl << "Quantity: " << quantity << endl;
        }
    private: //can only be accessed within class
        int calories;
}; //don't forget this semicolon!

//write this part later...(skip the first time through)
class turkey{
    public:
        void setTurkeyWeight(int x){ //mutator function - modifies the provate data member
            weight = x;
        }
        
        int getTurkeyWeight(){ // accessor function - just access the member without changing it
            return weight;
        }
        
        //skip these two lines the first time through...
        void setTurkeyName(string); //don't forget semicolon
        string getTurkeyName();
        //
    private:
        int weight;
        string name;
};
//skip until instructed
void turkey::setTurkeyName(string n){
    name = n;
}
string turkey::getTurkeyName(){
    return name;
}


int main(){
    
//classes with public and private members
    //classes = very similar to MATLAB structs.  Can group together variables of different types
    // C++ struct's are default public.  class defaults private.
        //Public class members - availiable outside of class
        //Private class members - can only be accessed only by the functions inside the class itself
        
//WRITE LINES 9-18 --------

    //declaration and initialization (just like any other variable type we want to use!)
    ThanksgivingFood food1; //declare the new class object
    food1.Name = "pumpkin pie"; //initialize the fields
    food1.quantity = 2;
    //run the displayInfo function from inside the class
    food1.displayInfo();
    //food1.calories = 100; // run this to show - will give compilation error. calories are private.
    
    //Typically, we want members to be private and we use a "public member function" to access it
//WRITE LINES 21-38
    turkey t1; //declare
    t1.setTurkeyWeight(12);
    cout << "This turkey weighs " << t1.getTurkeyWeight() << " pounds" << endl;
    
    //can also define those functions outside the class definition - "scope resolution operator ::"
//WRITE LINES 31-34 and 40-45
    t1.setTurkeyName("Gobble");
    cout << "This turkey's name is " << t1.getTurkeyName() << endl;
    
    
    //can also make arrays with the classes - array of objects
    turkey turk1[2];
    turk1[0].setTurkeyWeight(10);
    turk1[1].setTurkeyWeight(13);
    cout << "These turkeys weigh " << turk1[0].getTurkeyWeight() << " and " << turk1[1].getTurkeyWeight() << " pounds" << endl;
    
    
    
    // revisit how to compile C++ code in command line for the hw...
    // navigate to current folder with new terminal
    //type out:
    //        g++ -o run firstfile.cpp secondfile.cpp
    //hit enter.  in new line type:
    //        ./run
    //and hit enter again
    
    
    
    //if time: remind them how to read from csv file
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
