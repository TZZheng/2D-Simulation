import numpy as np
import matplotlib.pyplot as plt
import math
import random

#for Laplacian Solver
a = 2.8e-4
b = 5e-3
tau = .1
k = 0.005

size = 100
dx = 2./size

T = 1.0
dt = .9*dx**2/2
n = int(T/dt)

# U=np.random.rand(size,size)
# V=np.random.rand(size,size)
# def laplacian(Z):
#     Ztop=Z[0:-2,1:-1]
#     Zleft=Z[1:-1,0:-2]
#     Zbottom=Z[2:,1:-1]
#     Zright=Z[1:-1,2:]
#     Zcenter=Z[1:-1,1:-1]
#     return(Ztop+Zleft+Zbottom+Zright-4*Zcenter)/dx**2
# #for Mente Carlo Method
# for i in range(n):
#     deltaU= laplacian(U)
#     deltaV= laplacian(V)
#     Uc=U[1:-1,1:-1]
#     Vc=V[1:-1,1:-1]
#     U[1:-1,1:-1],V[1:-1,1:-1]=\
#         Uc+dt*(a*deltaU+Uc-Uc**3-Vc+k),\
#         Vc+dt*(b*deltaV+Uc-Vc)/tau
#     for Z in (U,V):
#         Z[0,:]=Z[1,:]
#         Z[-1,:]=Z[-2,:]
#         Z[:,0]=Z[:,1]
#         Z[:,-1]=Z[:,-2]
# plt.imshow(U,cmap=plt.cm.copper,extent=[-1,1,-1,1]);
# plt.xticks([]);
# plt.yticks([]);

#define parameters
E_field0 = 3.5e5
E_field = 3.5e5 #unit:V/cm
Eth_e = 1.2 #unit:eV
Eth_h = 1.5 #unit:eV
width_z = 1.6e-4 #unit: cm
width_x = 2e-2
width_y = 2e-2
v = 1e7 #assume maximum velocity,cm/s
q = 1.602e-19
t_resolution = 1e-14
D_z = 20 #cm^2/s
D_r = 10 #cm^2/s
R = 95000 #unit:Ohm/um^2

alpha_e = 1.286e6*math.exp(-1.4e6/E_field) # unit: /s POSSIBLE MODEL USED IN JIAN PAPER
alpha_h = 1.438e6*math.exp(-2.02e6/E_field)# unit: /s POSSIBLE MODEL USED IN JIAN PAPER
d_e = Eth_e/E_field
d_h = Eth_h/E_field
enabled_alpha_e = 1/(1/alpha_e - d_e)
enabled_alpha_h = 1/(1/alpha_h - d_h)
w_resolution = t_resolution*v # unit: cm
diff_step_z = math.sqrt(2*D_z*t_resolution) #units: cm % up to 6.3nm per step for diffusion up or down

def next_pos(Z, i): #unit:um
    Z[i, 2] = Z[i, 2]+(v*t_resolution+diff_step_z*random.random())*1e3
    Z[i, 0] = Z[i, 0]+(math.sqrt(2*D_r*t_resolution)*math.cos(random.random()*2*math.pi))*1e3
    Z[i, 1] = Z[i, 1]+(math.sqrt(D_r*t_resolution)*math.sin(random.random()*2*math.pi))*1e3
    return Z

total_simulation = 1e2
U_total = 5 #unit is not sure now
#rou = np.zeros(size, size, size)
jitter_dist = np.zeros((1, int(total_simulation)))
for i in range(int(total_simulation)):
    #set initial matrix
    pos_e = np.zeros((100000,3))
    t_e_next = np.zeros((100000, 1))
    if_active_e = np.zeros((100000, 1))
    pos_h = np.zeros((100000,3))
    t_h_next = np.zeros((100000, 1))
    t_h = np.zeros((100000, 1))
    if_active_h = np.zeros((100000, 1))
    total_carrier = 0
    counter = 0
    I = 0
    total_e = 1
    total_h = 1
    if_active_e[0] = 1
    if_active_h[0] = 1
    t_e_next[0]=d_e-math.log(random.random())/enabled_alpha_e/v
    # new U and rou matrix
    #U=PDE(U,rou)
    while(total_carrier< 1.25e4):
        alpha_e = 1.286e6*math.exp(-1.4e6/E_field) # unit: /s POSSIBLE MODEL USED IN JIAN PAPER
        alpha_h = 1.438e6*math.exp(-2.02e6/E_field) # unit: /s POSSIBLE MODEL USED IN JIAN PAPER
        d_e = Eth_e/E_field
        d_h = Eth_h/E_field
        enabled_alpha_e = 1/(1/alpha_e - d_e)
        enabled_alpha_h = 1/(1/alpha_h - d_h)
        counter = counter+1
        total_carrier_last = total_carrier
        for i in range(total_e):
            if (if_active_e[i] > 0):
                pos_e = next_pos(pos_e, i)
                if(pos_e[i][2]>width_z or abs(pos_e[i][0])>width_x*1e3 or abs(pos_e[i][0]>width_y*1e3)):
                    total_carrier=total_carrier-1
                    t_e_next[i]=-1
                    if_active_e[i]=0
                else:
                    if(t_e_next[i]<counter*t_resolution):
                        total_carrier=total_carrier+2
                        total_e=total_e+1
                        total_h=total_h+1
                        pos_e[total_e-1]=pos_e[i]
                        if_active_e[total_e-1]=1
                        pos_h[total_h-1]=pos_h[i]
                        if_active_h[total_h-1]=1
                        t_e_next[i]=((d_e-math.log(random.random())/enabled_alpha_e/v+counter*t_resolution))
                        t_e_next[total_e]=((d_e-math.log(random.random())/enabled_alpha_e/v+counter*t_resolution))
                        t_h_next[total_h]=((d_h-math.log(random.random())/enabled_alpha_h/v+counter*t_resolution))                               
        for i in range(total_h):
            if (if_active_h[i]!=0):
                pos_h=next_pos(pos_h,i)
                if(pos_h[i][2] > width_z or abs(pos_h[i][0])>width_x*1e3 or abs(pos_h[i][0]>width_y*1e3)):
                    total_carrier = total_carrier-1
                    t_h_next[i] = -1
                    if_active_h[i] = 0
                else:
                    if(t_h_next[i] > counter*t_resolution):
                        total_carrier = total_carrier+2
                        total_e = total_e+1
                        total_h=total_h+1
                        pos_e[total_e-1]=pos[i]
                        if_active_e[total_e-1]=1
                        pos_h[total_h-1]=pos[i]
                        if_active_h[total_h-1]=1
                        t_h_next[i]=(d_h-math.log(random())/enabled_alpha_h/v+counter*t_resolution)
                        t_e_next[total_e]=(d_e-math.log(random())/enabled_alpha_e/v+counter*t_resolution)
                        t_h_next[total_h]=(d_h-math.log(random())/enabled_alpha_h/v+counter*t_resolution)
        I=(total_carrier_last-total_carrier)*q/t_resolution
        pos_e_temp=pos_e[0:total_e]
        pos_h_temp=pos_h[0:total_h]
        max_e=pos_e_temp.max(axis=0)
        min_e=pos_e_temp.min(axis=0)
        max_h=pos_h_temp.max(axis=0)
        min_h=pos_h_temp.min(axis=0)
        max_x=max(max_e[0],max_h[0])
        min_x=min(min_e[0],min_h[0])
        max_y=max(max_e[1],max_h[1])
        min_y=min(min_e[1],min_h[1])
        square=(max_x-min_x)*(max_y-min_y)
        if(square!=0):
            E_field = (E_field0*width_z*1e-3-I*R/square)/(width_z*1e-3)
        if(counter%100==0):
            print(E_field)
            print(I," ",square)
            print(total_carrier)
    print(counter)
print(counter*t_resolution)
