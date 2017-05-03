#include<iostream>
#include<fstream>
#include<omp.h>
#include<math.h>
#include<time.h>
#include<malloc.h>
#define numarray 1000000
#define total_simulation 200
#define absorption_bin 100
#define random() (double)rand()/RAND_MAX
#define PI 3.1415926
#define max(a,b) (((a)>(b))?a:b)
#define min(a,b) (((a)<(b))?a:b)
#define threads 4
// another version using heap
using namespace std;
clock_t start_time;
clock_t end_time;
double max_array(double *Z, int i) {
	//find the max in an array
	double *p = Z;
	int count = 1;
	double max = *p;
	for (; count<i; ++count) {
		++p;
		if (*p>max) {
			max = *p;
		}
	}
	return max;
}
double min_array(double *Z, int i) {
	//find the min in an array
	double *p = Z;
	int count = 1;
	double min = *p;
	for (; count<i; ++count) {
		++p;
		if (*p<min) {
			min = *p;
		}
	}
	return min;
}
//double pos_e[numarray][3] = { 0 };
//double pos_h[numarray][3] = { 0 };
//int jitter_dist[total_simulation] = { 0 };
//double pos_e_temp0[numarray] = { 0 };
//double pos_e_temp1[numarray] = { 0 };
//double t_e_next[numarray] = { 0 };
//int if_active_e[numarray] = { 0 };
//double pos_h_temp0[numarray] = { 0 };
//double pos_h_temp1[numarray] = { 0 };
//double t_h_next[numarray] = { 0 };
//int  if_active_h[numarray] = { 0 };
int main() {
	//the physics parameters
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
	double D_r = 50; //cm^2/s
	double R = 75000; //unit:Ohm/um^2
	//the derivative theory parameters
	double alpha_e = 1286000 * exp(-1400000 / E_field); // unit: /s POSSIBLE MODEL USED IN JIAN PAPER
	double alpha_h = 1438000 * exp(-2020000 / E_field);// unit: /s POSSIBLE MODEL USED IN JIAN PAPER
	double d_e = Eth_e / E_field;
	double d_h = Eth_h / E_field;
	double enabled_alpha_e = 1 / (1 / alpha_e - d_e);
	double enabled_alpha_h = 1 / (1 / alpha_h - d_h);
	double w_resolution = t_resolution*v; // unit: cm
	double diff_step_z = sqrt(2 * D_z*t_resolution); //units: cm % up to 6.3nm per step for diffusion up or down
	double diff_step_r = sqrt(2 * D_r*t_resolution);

	//pointer for jitter_dist declaration
	int pointer = 0;
	// the absorption parameters
	double refraction_index = 3.65;
	double reflection_coef = pow(1 - refraction_index, 2) / pow(1 + refraction_index, 2);
	double absorption_coef = 535; //unit:cm
	double absorption[absorption_bin] = { 0 };

	// the loop-calculation definition
	//double total_carrier = 0;
	//int counter = 0;
	//double counter_I = 0;
	//double I = 0;
	//int total_e = 1;
	//int total_h = 1;
	int loop_temp = 0;
	int i, j, k, l, ii, iii, loop; //for loop
	double *array_pointer_e = NULL, *array_pointer_h = NULL;

	double min_x = 0;
	double max_x = 0;
	double min_y = 0;
	double max_y = 0;
	double square = 0;
	srand((unsigned int)time(NULL));
	// next use 12 parallel computing, 12 times space
	double *pos_e = NULL, *pos_h = NULL, *pos_e_temp0 = NULL, *pos_e_temp1 = NULL, *pos_h_temp0 = NULL, *pos_h_temp1 = NULL,*t_e_next=NULL,*t_h_next=NULL;
	//double *pos_e_pointer = NULL, *pos_h_pointer = NULL, *pos_e_temp0_pointer = NULL, *pos_e_temp1_pointer = NULL, *pos_h_temp0_pointer = NULL, *pos_h_temp1_pointer = NULL, *t_e_next_pointer = NULL, *t_h_next_pointer = NULL;
	int *if_active_e = NULL, *if_active_h = NULL;
	//int *if_active_e_pointer = NULL, *if_active_h_pointer = NULL;
	int counter_par = 0;
	pos_e = (double *)malloc(threads * 3 * numarray*sizeof(double));
	pos_h = (double *)malloc(threads * 3 * numarray*sizeof(double));
	if_active_e = (int *)malloc(threads * numarray*sizeof(int));
	if_active_h = (int *)malloc(threads * numarray*sizeof(int));
	int *jitter_dist = (int *)malloc(numarray*sizeof(int));
	pos_e_temp0 = (double *)malloc(threads * numarray*sizeof(double));
	pos_e_temp1 = (double *)malloc(threads * numarray*sizeof(double));
	pos_h_temp0 = (double *)malloc(threads * numarray*sizeof(double));
	pos_h_temp1 = (double *)malloc(threads * numarray*sizeof(double));
	t_e_next = (double *)malloc(threads * numarray*sizeof(double));
	t_h_next = (double *)malloc(threads * numarray*sizeof(double));
	for (loop = 0; loop<absorption_bin; ++loop) {
		absorption[loop] = (1 - reflection_coef)*(exp(-absorption_coef*(loop)*width_z / absorption_bin) - exp(-absorption_coef*(loop + 1)*width_z / absorption_bin)) + pow(1 - reflection_coef, 2)*(exp(-absorption_coef*(2 * width_z - (loop + 1)*width_z / absorption_bin)) - exp(-absorption_coef*(2 * width_z - loop*width_z / absorption_bin)));
		//cout<<"absorption"<<loop<<"\t"<<absorption[loop]<<endl;
	}
	for (i = 0; i<absorption_bin; ++i) {
		loop_temp = round(total_simulation*absorption[i]);
		counter_par = 0;
		//cout<<loop_temp<<endl;
		//memset(pos_e,0,sizeof(double)*3*numarray);
		//memset(t_e_next,0,sizeof(double)*numarray);
		//memset(if_active_e,0,sizeof(int)*numarray);
		//memset(pos_h,0,sizeof(double)*3*numarray);
		//memset(t_h_next,0,sizeof(double)*numarray);
		//memset(if_active_h,0,sizeof(int)*numarray);
		//start_time = clock();
		array_pointer_e = pos_e;
		array_pointer_h = pos_h;
		for (ii = 0; ii<numarray*3*threads; ++ii) {
			*array_pointer_e = 0;
			*array_pointer_h = 0;
			++array_pointer_e;
			++array_pointer_h;
		}
		array_pointer_e = t_e_next;
		array_pointer_h = t_h_next;
		for (ii = 0; ii<numarray*threads; ++ii) {
			++array_pointer_e;
			++array_pointer_h;
			if_active_e[ii] = 0;
			if_active_h[ii] = 0;
			t_e_next[ii] = 0;
			t_h_next[ii] = 0;
			pos_e_temp0[ii] = 0;
			pos_e_temp1[ii] = 0;
			pos_h_temp0[ii] = 0;
			pos_h_temp1[ii] = 0;
		}
		t_e_next[0] = (d_e - log(random())) / enabled_alpha_e / v;
		t_h_next[0] = (d_h - log(random())) / enabled_alpha_h / v;
		while (counter_par < loop_temp) {
			omp_set_num_threads(2);
			#pragma omp parallel for
			for (int loop_par = 0; loop_par < (loop_temp - counter_par>threads ? threads : loop_temp - counter_par); ++loop_par) {
				double *pos_e_pointer = pos_e + numarray * 3 * loop_par;
				double *pos_h_pointer = pos_h + numarray * 3 * loop_par;
				double *t_e_next_pointer = t_e_next + numarray * loop_par;
				double *t_h_next_pointer = t_h_next + numarray*loop_par;
				int *if_active_e_pointer = if_active_e + numarray*loop_par;
				int *if_active_h_pointer = if_active_h + numarray*loop_par;
				double *pos_e_temp0_pointer = pos_e_temp0 + numarray*loop_par;
				double *pos_e_temp1_pointer = pos_e_temp1 + numarray*loop_par;
				double *pos_h_temp0_pointer = pos_h_temp0 + numarray*loop_par;
				double *pos_h_temp1_pointer = pos_h_temp1 + numarray*loop_par;
				++counter_par;
				if_active_e_pointer[0] = 1;
				if_active_h_pointer[0] = 1;
				pos_e_pointer[0] = 0;
				pos_e_pointer[1] = 0;
				pos_e_pointer[2] = width_z / (absorption_bin)*i;
				pos_h_pointer[0] = 0;
				pos_h_pointer[1] = 0;
				pos_e_pointer[2] = width_z / (absorption_bin)*i;
				int total_carrier = 2;
				int total_e = 1;
				int total_h = 1;
				int counter = 0;
				int counter_I = 0;
				double I = 0;
				E_field = E_field0;
				while (total_carrier<12500 && total_carrier>0) {
					alpha_e = 1286000 * exp(-1400000 / E_field); // unit: /s POSSIBLE MODEL USED IN JIAN PAPER
					alpha_h = 1438000 * exp(-2020000 / E_field); // unit: /s POSSIBLE MODEL USED IN JIAN PAPER
					d_e = Eth_e / E_field;
					d_h = Eth_h / E_field;
					enabled_alpha_e = 1 / (1 / alpha_e - d_e);
					enabled_alpha_h = 1 / (1 / alpha_h - d_h);
					counter = counter + 1;
					j = 0;
					while(j<total_e) {
						if (if_active_e_pointer[j]>0) {
							pos_e_pointer[j * 3 + 2] = pos_e_pointer[j * 3 + 2] + (v*t_resolution + diff_step_z*random())*1e3;
							pos_e_pointer[j * 3] = pos_e_pointer[j * 3] + (diff_step_r*cos(random() * 2 * PI))*1e3;
							pos_e_pointer[j * 3 + 1] = pos_e_pointer[j * 3 + 1] + (diff_step_r*sin(random() * 2 * PI))*1e3;
							//TODO:  replace next function
							if (pos_e_pointer[j * 3 + 2]>width_z * 1000 || fabs(pos_e_pointer[j * 3])>width_x * 1000 || fabs(pos_e_pointer[j * 3 + 1]>width_y * 1000)) {
								--total_carrier;
								if (pos_e_pointer[j * 3 + 2] > width_z * 1000) {
									++counter_I;
								}
								t_e_next_pointer[j] = 999;
								if_active_e_pointer[j] = 0;
								pos_e_pointer[j * 3] = 0;
								pos_e_pointer[j * 3 + 1] = 0;
								pos_e_pointer[j * 3 + 2] = 0;
							}
							else {
								if (t_e_next_pointer[j]<counter*t_resolution) {
									total_carrier = total_carrier + 2;
									pos_e_pointer[total_e * 3] = pos_e_pointer[j * 3];
									pos_e_pointer[total_e * 3 + 1] = pos_e_pointer[j * 3 + 1];
									pos_e_pointer[total_e * 3 + 2] = pos_e_pointer[j * 3 + 2];
									if_active_e_pointer[total_e] = 1;
									pos_h_pointer[total_h * 3] = pos_h_pointer[j * 3];
									pos_h_pointer[total_h * 3 + 1] = pos_h_pointer[j * 3 + 1];
									pos_h_pointer[total_h * 3 + 2] = pos_h_pointer[j * 3 + 2];
									if_active_h_pointer[total_h] = 1;
									t_e_next_pointer[j] = (d_e - log(random())) / enabled_alpha_e / v + counter*t_resolution;
									t_e_next_pointer[total_e] = (d_e - log(random())) / enabled_alpha_e / v + counter*t_resolution;
									t_h_next_pointer[total_h] = (d_h - log(random())) / enabled_alpha_h / v + counter*t_resolution;
									++total_e;
									++total_h;
								}
							}
						}
						++j;
					}
					j = 0;
					while(j<total_h) {
						if (if_active_h_pointer[j]>0) {
							pos_h_pointer[j * 3 + 2] = pos_h_pointer[j * 3 + 2] - (v*t_resolution + diff_step_z*random())*1e3;
							pos_h_pointer[j * 3] = pos_h_pointer[j * 3] + (sqrt(2 * D_r*t_resolution)*cos(random() * 2 * PI))*1e3;
							pos_h_pointer[j * 3 + 1] = pos_h_pointer[j * 3 + 1] + (sqrt(2 * D_r*t_resolution)*sin(random() * 2 * PI))*1e3;
							//TODO:  replace next function
							if (pos_h_pointer[j * 3 + 2]<0 || fabs(pos_h_pointer[j * 3])>width_x * 1000 || fabs(pos_h_pointer[j * 3 + 1]>width_y * 1000)) {
								total_carrier = total_carrier - 1;
								if (pos_h_pointer[j * 3 + 2] > width_z * 1000) {
									counter_I = counter_I + 1;
								}
								t_h_next_pointer[j] = 999;
								if_active_h_pointer[j] = 0;
								pos_h_pointer[j * 3] = 0;
								pos_h_pointer[j * 3 + 1] = 0;
								pos_h_pointer[j * 3 + 2] = 0;
							}
							else {
								if (t_h_next_pointer[j]<counter*t_resolution) {
									total_carrier = total_carrier + 2;
									pos_e_pointer[total_e * 3] = pos_e_pointer[j * 3];
									pos_e_pointer[total_e * 3 + 1] = pos_e_pointer[j * 3 + 1];
									pos_e_pointer[total_e * 3 + 2] = pos_e_pointer[j * 3 + 2];
									if_active_e_pointer[total_e] = 1;
									pos_h_pointer[total_h * 3] = pos_h_pointer[j * 3];
									pos_h_pointer[total_h * 3 + 1] = pos_h_pointer[j * 3 + 1];
									pos_h_pointer[total_h * 3 + 2] = pos_h_pointer[j * 3 + 2];
									if_active_h_pointer[total_h] = 1;
									t_h_next_pointer[j] = (d_h - log(random())) / enabled_alpha_h / v + counter*t_resolution;
									t_e_next_pointer[total_e] = (d_e - log(random())) / enabled_alpha_e / v + counter*t_resolution;
									t_h_next_pointer[total_h] = (d_h - log(random())) / enabled_alpha_h / v + counter*t_resolution;
									total_e = total_e + 1;
									total_h = total_h + 1;
								}
							}
						}
						++j;
					}
					I = counter_I*q / t_resolution;
					//find the max and min value to calculate square
					if (total_carrier>200) {
						k = 0;
						while( k<total_e) {
							pos_e_temp0_pointer[k] = pos_e_pointer[k * 3];
							pos_e_temp1_pointer[k] = pos_e_pointer[k * 3 + 1];
							++k;
						}
						k = 0;
						while( k<total_h) {
							pos_h_temp0_pointer[k] = pos_h_pointer[k * 3];
							pos_h_temp1_pointer[k] = pos_h_pointer[k * 3 + 1];
							++k;
						}
						max_x = max(max_array(pos_e_temp0_pointer, total_e), max_array(pos_h_temp0_pointer, total_h));
						max_y = max(max_array(pos_e_temp1_pointer, total_e), max_array(pos_h_temp1_pointer, total_h));
						min_x = min(min_array(pos_e_temp0_pointer, total_e), min_array(pos_h_temp0_pointer, total_h));
						min_x = min(min_array(pos_e_temp1_pointer, total_e), min_array(pos_h_temp1_pointer, total_h));
						square = (max_x - min_x)*(max_y - min_y);
						E_field = E_field0 - I*R / square / (width_z * 1000);
						//cout<<I<<"\t"<<square<<"\t"<<E_field<<endl;
						cout << omp_get_thread_num()<< endl;
					}
				}
				if (total_carrier > 0) {
					jitter_dist[pointer] = counter;
				}
				else {
					jitter_dist[pointer] = 0;
				}
				cout << total_carrier << endl;
				//end_time = clock();
				//cout << "Running time is: " << static_cast<double>(end_time - start_time) / CLOCKS_PER_SEC * 1000 << "ms" << endl;
				pointer = pointer + 1;
			}
		}
		cout << i << endl;
	}
	ofstream SaveFile("prototype.txt");
	//int temp;
	for (iii = 0; iii<pointer; ++iii) {
		SaveFile << jitter_dist[iii] << endl;
	}
	SaveFile.close();
	return 0;
}