#include <iostream>
#include <cstdlib>
#include <cmath>

using namespace std;

int main()
{
    
    // Variables declaration and Initialization//	   
    const double m  = 2 * 1.00784 * 1.660538921e-27; // H2 unit mass [kg]
    const double kb = 1.38064852e-23;                // Boltzmann’s constant [J/K]
    const double na = 6.02214076e+23;                // Avogadro's constant
    const double c  = 1000.0;                        // Specific heat capacity [J/kg/K]
    int n;
    double vmin, vmax, Ts, Tg, Q;
   
    // Ask for data entry from console
    cout << "Enter the number of particles: ";
    cin  >> n;
    cout << "Enter the minimum velocity in m/s: ";
    cin  >> vmin;
    cout << "Enter the maximum velocity in m/s: ";
    cin  >> vmax;
    cout << "Enter the temperature of surrounding air in K: ";
    cin  >> Ts;

    // Compute Tg (Gas Temperature) by iteration - for-loop from i=1 to n
    srand (1);                                // initialize the random number generator seed 
    double vsqsum = 0;                        // initialize sum of the squared velocity
    double velocity, r;                       // declare loop-body variables 
    for (int i = 1; i <= n; i++)       
    {
       r = rand() * (1.0 / RAND_MAX);         // generate random number between 0 and 1
       velocity = r * (vmax - vmin) + vmin;   // scale r to values between vmin and vmax
       vsqsum = vsqsum + velocity * velocity;
    }
    
    // compute and output the gas temperature to console
    Tg = m * vsqsum / n / (3 * kb);          // Eq (1) 
    cout << "The gas temperature in K: " << round(Tg) << endl;
    
    // Compute heat transfer amount 
    Q = round(c * na * m * (Ts - Tg));              // Eq (2) 
    
    // Output respective messages according to Ts - Tg conditions
    if (Q == 0) 
    {
        cout << "No heat transfer" << endl;
    }
    else if (Q > 0) 
    {
        cout << "Heat is transferred from the surrounding air to the gas" << endl;
        cout << "The transferred heat in J: " << Q << endl;
    }
    else 
    {
        cout << "Heat is transferred away from the gas to the surrounding air" << endl;
        cout << "The transferred heat in J: " << Q << endl;
    }  
    
    return 0; // return with no error
}