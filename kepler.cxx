#include <iostream>
#include <fstream> // include this for output file
#include <math.h> // include this for pi, pow() and ceil()

using namespace std;


void ellips(double* p, double* q, const double dt); // function to calc momentum and position

int main(){
    
    double p[2]; //declare momentum
    double q[2]; // declare position
    double H; // declase Hamiltonian
    
    
    double const t_end = 20*M_PI; // max. simulation time
    double const dt = 0.0005; // time step
    const int N = ceil(t_end/dt); // max. amount of steps
    double const e = 0.6; // some constant
    
    const string st = to_string(dt); // string to rename output files
    const string filename = "data_kepler" "_" + st + ".txt"; // make output file name
    ofstream out(filename.c_str());
    
    
    p[0] = 0  ; p[1] = sqrt((1.0 + e) / (1.0 - e)) ; q[0] = 1.0 - e ; q[1] = 0 ; // initial conditions
    
    
    H = 0.5 * (pow( p[0],2) + pow(p[1],2)) - 1.0 / sqrt(pow(q[0],2) + pow(q[1],2)); // Hamiltonian
    
    
    out << 0 << "\t" << H << "\t" << q[0] << "\t" << q[1] << endl; // write data for t = 0

    
    for (int i = 1; i< N-1; i++) { // calculate momentua and positions
        ellips(p,q,dt);
        
        H = 0.5 * (pow( p[0],2) + pow(p[1],2)) - 1.0 / sqrt(pow(q[0],2) + pow(q[1],2));
    
        out << i*dt << "\t" << H << "\t" << q[0] << "\t" << q[1] << endl; // write all values to output file
    
    
    }
    out.close(); // close output file
    return 0;
}


void ellips(double* p, double* q, const double dt){ //function to calc momentum and position

    p[0] = p[0] - dt * q[0] / pow(pow(q[0],2) + pow(q[1],2),1.5); // momenta
    p[1] = p[1] - dt * q[1] / pow(pow(q[0],2) + pow(q[1],2),1.5);
    
    q[0] = q[0] + dt * p[0]; // positions
    q[1] = q[1] + dt * p[1];

}