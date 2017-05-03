#include<iostream>
#include<math.h>
#include<time.h>
#define numarray 1000000
#define total_simulation 10000
#define absorption_bin 100
#define random() (double)rand()/RAND_MAX
#define PI 3.1415926
using namespace std;

double max_array(double *Z, int i){
	double *p=Z;
	int count=1;
	double max=*p;
	for(;count<i;++count){
		++p;
		if(*p>max){
			max=*p;
		}
	}
	return max;
}
double min_array(double *Z, int i){
	double *p=Z;
	int count=1;
	double min=*p;
	for(;count<i;++count){
		++p;
		if(*p<min){
			min=*p;
		}
	}
	return min;
}    
double pos_e[numarray][3]={0};
double pos_h[numarray][3]={0};
int jitter_dist[total_simulation]={0};
double pos_e_temp0[numarray]={0};
double pos_e_temp1[numarray]={0};
double t_e_next[numarray]={0};
int if_active_e[numarray]={0};
double pos_h_temp0[numarray]={0};
double pos_h_temp1[numarray]={0};
double t_h_next[numarray]={0};
int  if_active_h[numarray]={0};
int main(){
    double E_field0 = 3.5e5;
    double E_field = 3.5e5; //unit:V/cm
    double Eth_e = 1.2; //unit:eV
    double Eth_h = 1.5; //unit:eV
    double width_z = 1.6e-3; //unit: cm
    double width_x = 2e-002;
    double width_y = 2e-002;
    double v = 1e7; //assume maximum velocity,cm/s
    double q = 1.602e-019;
    double t_resolution = 1e-014;
    double D_z = 20; //cm^2/s
    double D_r = 0; //cm^2/s
    double R = 0; //unit:Ohm/um^2

    double alpha_e = 1.286e6*exp(-1.4e6/E_field); // unit: /s POSSIBLE MODEL USED IN JIAN PAPER
    double alpha_h = 1.438e6*exp(-2.02e6/E_field);// unit: /s POSSIBLE MODEL USED IN JIAN PAPER
    double d_e = Eth_e/E_field;
    double d_h = Eth_h/E_field;
    double enabled_alpha_e = 1/(1/alpha_e - d_e);
    double enabled_alpha_h = 1/(1/alpha_h - d_h);
    double w_resolution = t_resolution*v; // unit: cm
    double diff_step_z = sqrt(2*D_z*t_resolution); //units: cm % up to 6.3nm per step for diffusion up or down


    int pointer=0;

    double refraction_index=3.65;
    double reflection_coef=pow(1-refraction_index,2)/pow(1+refraction_index,2);
    double absorption_coef=535; //unit:cm
    double absorption[absorption_bin]={0};


    double total_carrier=0;
    double counter = 0;
    double counter_I = 0;
    double I = 0;
    int total_e = 1;
    int total_e_temp=1;
    int total_h = 1;
    int total_h_temp=1;
    int loop_temp=0;
    
	double min_x=0;
	double max_x=0;
	double min_y=0;
	double max_y=0;
	double square=0;
    srand((unsigned int)time(NULL));
    for(int loop=0;loop<absorption_bin;++loop){
        absorption[loop]=(1-reflection_coef)*(exp(-absorption_coef*(loop)*width_z/absorption_bin)-exp(-absorption_coef*(loop+1)*width_z/absorption_bin))+pow(1-reflection_coef,2)*(exp(-absorption_coef*(2*width_z-(loop+1)*width_z/absorption_bin))-exp(-absorption_coef*(2*width_z-loop*width_z/absorption_bin)));
        //cout<<"absorption"<<loop<<"\t"<<absorption[loop]<<endl;
    }
    for(int i=0;i<absorption_bin;++i){
		loop_temp=round(total_simulation*absorption[i]);
        //cout<<loop_temp<<endl;
        for(int loop=0;loop<loop_temp;++loop){
            //memset(pos_e,0,sizeof(double)*3*numarray);
            //memset(t_e_next,0,sizeof(double)*numarray);
            //memset(if_active_e,0,sizeof(int)*numarray);
            //memset(pos_h,0,sizeof(double)*3*numarray);
            //memset(t_h_next,0,sizeof(double)*numarray);
            //memset(if_active_h,0,sizeof(int)*numarray);
            for(int ii=0;ii<numarray;++ii){
            	pos_e[ii][0]=0;
            	pos_e[ii][1]=0;
            	pos_e[ii][2]=0;
            	pos_h[ii][0]=0;
            	pos_h[ii][1]=0;
            	pos_h[ii][2]=0;
            }
            for(int ii=0;ii<numarray;++ii){
            	t_e_next[ii]=0;
            	t_h_next[ii]=0;
            	if_active_e[ii]=0;
            	if_active_e[ii]=0;
            }
            if_active_e[0] = 1;
            if_active_h[0] = 1;
            t_e_next[0] = (d_e-log(random()))/enabled_alpha_e/v;
            t_h_next[0] = (d_h-log(random()))/enabled_alpha_h/v;
            pos_e[0][0] = 0;
            pos_e[0][1] = 0;
            pos_e[0][2] = width_z/(absorption_bin)*i;
            pos_h[0][0] = 0;
            pos_h[0][1] = 0;
            pos_h[0][2] = width_z/(absorption_bin)*i;
            total_carrier=2;
            total_e=1;
			total_h=1; 
            counter=0;
            counter_I=0;
            E_field = E_field0;
            while(total_carrier<12500 && total_carrier>0){
                alpha_e = 1286000*exp(-1400000/E_field); // unit: /s POSSIBLE MODEL USED IN JIAN PAPER
                alpha_h = 1438000*exp(-2020000/E_field); // unit: /s POSSIBLE MODEL USED IN JIAN PAPER
                d_e = Eth_e/E_field;
                d_h = Eth_h/E_field;
                enabled_alpha_e = 1/(1/alpha_e - d_e);
                enabled_alpha_h = 1/(1/alpha_h - d_h);
                counter = counter+1;
                total_e_temp=total_e;
                for(int j=0;j<total_e_temp;++j){
                    if(if_active_e[j]>0){
                        pos_e[j][2] += (v*t_resolution+diff_step_z*(2*random()-1));
                        //pos_e[j][0] = pos_e[j][0]+(sqrt(2*D_r*t_resolution)*cos(random()*2*PI))*1e3;
                        //pos_e[j][1] = pos_e[j][1]+(sqrt(2*D_r*t_resolution)*sin(random()*2*PI))*1e3;
                        //TODO:  replace next function
                        if(pos_e[j][2]>width_z or fabs(pos_e[j][0])>width_x or fabs(pos_e[j][1]>width_y)){
                            --total_carrier;
                            if (pos_e[j][2] > width_z){
                                ++counter_I;
                            }
                            if(fabs(pos_e[j][0])>width_x || fabs(pos_e[j][1]>width_y)){
                            	cout<<"bang!"<<endl;
                            }
                            t_e_next[j] = 999;
                            if_active_e[j] = 0;
                            pos_e[j][0] = 0;
                            pos_e[j][1] = 0;
                            pos_e[j][2] = 0;
                        }
                        else{
                            if(t_e_next[j]<counter*t_resolution){
                                total_carrier=total_carrier+2;
                                pos_e[total_e][0]=pos_e[j][0];
                                pos_e[total_e][1]=pos_e[j][1];
                                pos_e[total_e][2]=pos_e[j][2];
                                if_active_e[total_e]=1;
                                pos_h[total_h][0]=pos_h[j][0];
                                pos_h[total_h][1]=pos_h[j][1];
                                pos_h[total_h][2]=pos_h[j][2];
                                if_active_h[total_h]=1;
                                t_e_next[j]=(d_e-log(random()))/enabled_alpha_e/v+counter*t_resolution;
                                t_e_next[total_e]=(d_e-log(random()))/enabled_alpha_e/v+counter*t_resolution;
                                t_h_next[total_h]=(d_h-log(random()))/enabled_alpha_h/v+counter*t_resolution;
                                ++total_e;
                                ++total_h;
                            }
                        }
                    }
                }
                total_h_temp=total_h;
                for(int j=0;j<total_h_temp;++j){
                    if(if_active_h[j]>0){
                        pos_h[j][2] = pos_h[j][2]-(v*t_resolution+diff_step_z*(2*random()-1));
                        pos_h[j][0] = pos_h[j][0]+(sqrt(2*D_r*t_resolution)*cos(random()*2*PI));
                        pos_h[j][1] = pos_h[j][1]+(sqrt(2*D_r*t_resolution)*sin(random()*2*PI));
                        //TODO:  replace next function
                        if(pos_h[j][2]<0 || fabs(pos_h[j][0])>width_x || fabs(pos_h[j][1]>width_y)){
                            total_carrier=total_carrier-1;
                            if (pos_h[j][2] > width_z){
                                counter_I = counter_I+1;
                            }
                            if(fabs(pos_h[j][0])>width_x || fabs(pos_h[j][1]>width_y)){
                            	cout<<"bang!"<<endl;
                            }
                            t_h_next[j] = 999;
                            if_active_h[j] = 0;
                            pos_h[j][0] = 0;
                            pos_h[j][1] = 0;
                            pos_h[j][2] = 0;
                        }
                        else{
                            if(t_h_next[j]<counter*t_resolution){
                                total_carrier=total_carrier+2;
                                pos_e[total_e][0]=pos_e[j][0];
                                pos_e[total_e][1]=pos_e[j][1];
                                pos_e[total_e][2]=pos_e[j][2];
                                if_active_e[total_e]=1;
                                pos_h[total_h][0]=pos_h[j][0];
                                pos_h[total_h][1]=pos_h[j][1];
                                pos_h[total_h][2]=pos_h[j][2];
                                if_active_h[total_h]=1;
                                t_h_next[j]=(d_h-log(random()))/enabled_alpha_h/v+counter*t_resolution;
                                t_e_next[total_e]=(d_e-log(random()))/enabled_alpha_e/v+counter*t_resolution;
                                t_h_next[total_h]=(d_h-log(random()))/enabled_alpha_h/v+counter*t_resolution;
                                total_e=total_e+1;
                                total_h=total_h+1;
                            }
                        }
                    }
                }
                I=counter_I*q/t_resolution;
                //find the max and min value to calculate square
                if(total_carrier>20000){
                	for(int k=0;k<total_e;k++){
                    	pos_e_temp0[k]=pos_e[k][0];
                    	pos_e_temp1[k]=pos_e[k][1];
                	}
                	for(int k=0;k<total_h;k++){
                    	pos_h_temp0[k]=pos_h[k][0];
                    	pos_h_temp1[k]=pos_h[k][1];
                	}
                	max_x=max(max_array(pos_e_temp0,total_e),max_array(pos_h_temp0,total_h));
                	max_y=max(max_array(pos_e_temp1,total_e),max_array(pos_h_temp1,total_h));
                	min_x=min(min_array(pos_e_temp0,total_e),min_array(pos_h_temp0,total_h));
                	min_x=min(min_array(pos_e_temp1,total_e),min_array(pos_h_temp1,total_h));
                	square=(max_x-min_x)*(max_y-min_y);
                    E_field = E_field0-I*R/square/(width_z);
                    //cout<<I<<"\t"<<square<<"\t"<<E_field<<endl;
                }
            }
            if (total_carrier>12500){
                jitter_dist[pointer] = counter;
            }
            else{
            	jitter_dist[pointer] = 0;
            } 
            pointer=pointer+1;
            //cout<<counter<<endl;
        }
        cout<<i<<endl;
    }
    FILE *fp=fopen("prototype.txt","w");
    int temp;
    for(int iii=0;iii<pointer;++iii){
    	fprintf(fp,"%d\n",jitter_dist[iii]);	
    }
    fclose(fp);
    return 0;
}
