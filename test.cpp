#include<iostream>
#include<math.h>
#include<memory.h>
#include<time.h>
#define numarray 100000
#define total_simulation 10000
#define absorption_bin 100
#define random() (double)rand()/RAND_MAX
#define PI 3.1415926
using namespace std;    
double pos_e[numarray][3]={0};
double pos_h[numarray][3]={0};
int main(){
	double E_field0 = 350000;
    double E_field = 350000; //unit:V/cm
    double Eth_e = 1.2; //unit:eV
    double Eth_h = 1.5; //unit:eV
    double width_z = 1.6e-003; //unit: cm
    double width_x = 2e-002;
    double width_y = 2e-002;
    double v = 10000000; //assume maximum velocity,cm/s
    double q = 1.602e-019;
    double t_resolution = 1e-014;
    double D_z = 20; //cm^2/s
    double D_r = 100; //cm^2/s
    double R = 75000; //unit:Ohm/um^2

    double alpha_e = 1286000*exp(-1400000/E_field); // unit: /s POSSIBLE MODEL USED IN JIAN PAPER
    double alpha_h = 1438000*exp(-2020000/E_field);// unit: /s POSSIBLE MODEL USED IN JIAN PAPER
    double d_e = Eth_e/E_field;
    double d_h = Eth_h/E_field;
    double enabled_alpha_e = 1/(1/alpha_e - d_e);
    double enabled_alpha_h = 1/(1/alpha_h - d_h);
    double w_resolution = t_resolution*v; // unit: cm
    double diff_step_z = sqrt(2*D_z*t_resolution); //units: cm % up to 6.3nm per step for diffusion up or down

    int jitter_dist[total_simulation]={0};
    int pointer=0;

    double refraction_index=3.65;
    double reflection_coef=pow(1-refraction_index,2)/pow(1+refraction_index,2);
    double absorption_coef=535; //unit:cm
    double absorption[absorption_bin]={0};
    
    
    double pos_e_temp0[numarray];
    double pos_e_temp1[numarray];
    double t_e_next[numarray];
    int if_active_e[numarray];
    double pos_h_temp0[numarray];
    double pos_h_temp1[numarray];
    double t_h_next[numarray];
    int  if_active_h[numarray];
    double total_carrier=0;
    double counter = 0;
    double counter_I = 0;
    double I = 0;
    int total_e = 1;
    int total_h = 1;
    
	double min_x=0;
	double max_x=0;
	double min_y=0;
	double max_y=0;
	double square=0;
	cout<<I;
}
