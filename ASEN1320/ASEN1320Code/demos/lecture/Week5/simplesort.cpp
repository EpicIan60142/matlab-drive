#include<iostream>

using namespace std;

int main (){
  
   int temp;
   int a[10] = {10,2,0,14,43,25,18,1,5,45};
   
   cout << "Before sorting: ";
   for(int i=0; i<10; i++) cout<<" "<<a[i];
   cout << endl;
  
   for(int i = 0; i<10; i++){
   for(int j = i+1; j<10; j++){
      if(a[j] < a[i]) 
      {
         temp = a[i];
         a[i] = a[j];
         a[j] = temp;
      }
   }}

   cout <<"After sorting:  ";
   for(int i=0; i<10; i++) cout<<" "<<a[i];
   cout << endl; 
   
   return 0;
}