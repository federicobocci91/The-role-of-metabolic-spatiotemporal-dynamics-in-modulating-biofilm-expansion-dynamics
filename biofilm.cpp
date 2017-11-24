#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <time.h>
#include <math.h>

using namespace std;


// function to compute first derivative
void der(int size, double dr, vector< double > &Z, vector< double > &der_Z)
{
    for(int i = 1; i<size-1; i++)
    {
        der_Z[i]= (Z[i+1]-Z[i-1])/(2*dr);
    }
}


// function to compute second derivative
void sec_der(int size, double dr, vector<double> &Z, vector<double> &lap_Z)
{
    for (int i = 1; i<size-1; i++) {
        lap_Z[i-1]= (Z[i+1]+Z[i-1]-2.*Z[i])/(dr*dr);
    }
}


// boundary conditions (at boundary Glutamate = G_ext and anything else is zero)
void dir_BC(int size, vector<double> &A,vector<double> &G,vector<double> &H, vector<double> &R, double G_ext, double g_0)
{
	// continuity condition on the first point of radial grid
    G[0]=G[1];
    A[0]=A[1];
    H[0]=H[1];
    R[0]=R[1];
	
	// biofilm boundary
    G[size-1]=G_ext/g_0;
    A[size-1]=0.;
    H[size-1]=0.;
    R[size-1]=0.;
    
}


// this function computes the sigmoid signal (1D radial distribution)
void sigmoid(vector<double> &R,vector<double> &f, int size, double r_s, double l_0)
{
    for (int i  = 0; i<size; i++) {
        f[i] = 1. - 1./(1. + exp(-(R[i]-r_s)/l_0));
    }
}


// compute amount of housekeeping proteins for biofilm growth
double R_growth(vector<double> &R, int size, double dr, double R0)
{
    double n = 0.;
    for (int i = 0; i<size; i++) {
        if (R[i]-R0>0) {
            n = n + 6.28*(dr)*(i+1)*(dr)*(R[i]-R0);
        }
    }
    return n;
}




// program starts here

int main()
{
	
	double G_ext = 30.;           		// external glutamate concentration   (mM)
	double m =12.;                      // hill coefficient for GDH activation
	double K_h = 3.;                    // glutamate threshold for GDH activation  (mM) // 1
	double D_A =  10000.;               // diffusion coefficient for Ammonia   (microm^2 h^-1)
	double D_G = 0.1*D_A; //            // diffusion constant for Glutamate  (microm^2 h^-1)
	double b_h = 50.*pow(10.,-1.);      // GHD activation coefficient  (mM h^-1)
	double a = 50.*pow(10.,1.);         // Ammonia production rate  (mM^-1 h^-1)
	double d_a = 4.*pow(10.,3.) ;       // consumption rate of Ammonia from R  (mM^-1 h^-1)
	double d_g = 3.2*pow(10.,3.);		// consumption rate of Glutamate from R  (mM^-1 h^-1)
	double b_r = 0.14;                  // Housekeeping proteins production rate  (mM^-1 h^-1)
	double g_h = 8.;                    // GDH degradation rate  (h^-1)
	double g_r = 1.6; // 1.6            // Housekeeping proteins degradation rate (h^-1)
	
	double f = pow(K_h,m);
	
	// parameter of dimensionless model
	double D = D_G/D_A;
	double alfa = a*b_h/(g_h*g_h);
	double delta = d_g/d_a;
	double beta = a*K_h*K_h*b_h*b_r*d_a/pow(g_h,4.);
	double g = g_r/g_h;
	
	// scaling factors for dimensionless model
	double tau = 1./g_h;
	double l = sqrt(D_A/g_h);
	double a_0 = a*K_h*b_h/(g_h*g_h);
	double g_0 = K_h;
	double h_0 = b_h/g_h;
	double r_0 = g_h/d_a;
	
	double b_height = 6.;             // biofilm height
    double k = 11.25*b_height/D_A;   // k = 11.25; 6 = biofilm height
    double R0=0.0001/r_0;             // housekeeping proteins threshold for growth

	// shaprness of interior-periphery transition
	double l_0 = 1.5/l;
	
	// integrator set-up
	int size = 100;
	double dr =  1.0/l;
	double T = 1000.;
	double dt =  0.00025;
	int n = int(T/dt);

	// peripheral layer parameters
	double G_threshold = 15.;       // glutamate threshold for peripheral layer
	double d0 = 150./l;             // saturation constant
	double const_d = d0*( (G_ext/g_0)/( (G_ext/g_0)+ (G_threshold/g_0)) );
	
	double const_s = 450.;
	double d_0 = const_d*tanh(size*dr*l/(const_s));  // width of peripheral layer
	double r_s = size*dr - d_0;                      // radius of interior region

	// vectors for all quantities and their first and second derivative
	vector< double > A(size);
	vector< double > G(size);
	vector< double > H(size);
	vector< double > R(size);
	vector< double > Radius(size);
	vector< double > sigm(size);
    
	for (int y =0; y<size; y++) {
		A[y]=0.;
		G[y]=0.;
		H[y]=0.;
		R[y]=0.;
		Radius[y] = dr*(y+0.5);   // vector of radial coordinate
	}
	sigmoid(Radius,sigm, size,r_s,l_0);

	vector< double > A_sec_der(size-2);
	vector< double > G_sec_der(size-2);
	vector< double > A_der(size-2);
	vector< double > G_der(size-2);

	ofstream myfile;
	myfile.open("data.txt");
    
	for (int i = 0; i<n; i++) {
		
		// boundary conditions
		dir_BC(size,A,G,H,R, G_ext, g_0);
			
		// compute first and second derivative for ammonium (A) and glutamate (G):
		sec_der(size,dr,A,A_sec_der);
		sec_der(size,dr,G,G_sec_der);
		der(size,dr,A,A_der);
		der(size,dr,G,G_der);
			
		////////////////////////////////////
		//                                //
		//    metabolic network eqs.      //
		//                                //
		////////////////////////////////////
        
		for (int j = 0; j<size-2; j++) {
            
			A[j+1] = A[j+1] + dt*( (A_sec_der[j]+ (A_der[j]/Radius[j+1])) + G[j+1]*H[j+1] - A[j+1]*R[j+1]);
			G[j+1] = G[j+1] + dt*( D*(G_sec_der[j]+ (G_der[j]/Radius[j+1])) - alfa*G[j+1]*H[j+1] - delta*G[j+1]*R[j+1]);
			H[j+1] = H[j+1] + dt*( (pow(G[j+1],m)/(pow(G[j+1],m) + 1.))*sigm[j+1] - H[j+1]);
			R[j+1] = R[j+1] + dt*( beta*G[j+1]*A[j+1] - g*R[j+1] );
		}
    
		////////////////////////////////////
		//                                //
		//       growth algorithm         //
		//                                //
		////////////////////////////////////
		
		if ( i*dt>112.5) {    // initial equilibration time (112.5)
				 
		// compute increment of radius
			double inc = ((k)*dt*(R_growth(R,size,dr,R0)))/(size*size*dr);
			dr = dr + inc;
			double scaling_factor = pow((dr-inc)/dr,2.);  // scaling factor (old area)/(new area)
				 
			// update radial coordinate vector and all variables
			for (int y =0; y<size; y++) {
				Radius[y] = dr*(y+0.5);
				A[y] = A[y]*scaling_factor;
				G[y] = G[y]*scaling_factor;
				H[y] = H[y]*scaling_factor;
				R[y] = R[y]*scaling_factor;
			}
				 
			// compute new peripheral layer, interior radius and update sigmaid signal
			d_0 = const_d*tanh(size*dr*l/(const_s));
			r_s = size*dr - d_0;
			sigmoid(Radius,sigm, size,r_s,l_0);
				 
			// save time, radius and growth rate to file data.txt
			if (i%100==0) {
				myfile << (i*dt-112.5)*tau << "   " << dr*l*size << "   " << l*(inc*size)/(dt*tau) << endl;
			}
		}
	}
    
	dir_BC(size,A,G,H,R, G_ext, g_0);
		
	myfile.close();

}
