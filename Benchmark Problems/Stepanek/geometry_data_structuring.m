function [X,Y,dx,dy,sigma_t,sigma_s,nu_sigma_f]=geometry_data_structuring()

%region 1
sigma_t_1=0.6;
sigma_s_1=0.53;
nu_sigma_f_1=0.079;

%region 2

sigma_t_2=0.48;
sigma_s_2=0.20;
nu_sigma_f_2=0;

%region 3
sigma_t_3=0.70;
sigma_s_3=0.66;
nu_sigma_f_3=0.043;

%region 4
sigma_t_4=0.65;
sigma_s_4=0.50;
nu_sigma_f_4=0;

%region 5
sigma_t_5=0.90;
sigma_s_5=0.89;
nu_sigma_f_5=0;

%% geometry

X=96;
Y=86;

dx=0.1;
dy=0.1;

x=(0:dx:X)';
y=(0:dy:Y)';

mesh_count_x=length(x)-1;
mesh_count_y=length(y)-1;

%% initialization
sigma_t=zeros(mesh_count_y,mesh_count_x);
sigma_s=zeros(mesh_count_y,mesh_count_x);
nu_sigma_f=zeros(mesh_count_y,mesh_count_x);

%%

% region 1
x_s_1=18;
y_s_1=18;

x_e_1=x_s_1+30;
y_e_1=y_s_1+25;

% region 2
x_s_2=x_e_1;;
y_s_2=18;

x_e_2=x_s_2+30;
y_e_2=y_s_2+25;

% region 3
x_s_3=x_e_1;;
y_s_3=18+y_e_2;

x_e_3=x_s_3+30;
y_e_3=y_s_3+25;

% region 4
x_s_4=18;;
y_s_4=18+y_e_1;

x_e_4=x_s_4+30;
y_e_4=y_s_4+25;

% region 5

x_s_5=0;
y_s_5=0;

x_e_5=96;
y_e_5=86;

%% mesh geometry

% region 5 section 1 and 2
sec_1_x=(0:dx:X)';
sec_1_y=(0:dy:18)';

mesh_p_sec_1_x=length(sec_1_x)-1;
mesh_p_sec_1_y=length(sec_1_y)-1;

sec_2_x=(0:dx:18)';
sec_2_y=(18:dy:Y-18)';

mesh_p_sec_2_x=length(sec_2_x)-1;
mesh_p_sec_2_y=length(sec_2_y)-1+mesh_p_sec_1_y;

%region 1

r_1_x=(x_s_1:dx:x_e_1);
r_1_y=(y_s_1:dy:y_e_1);

mesh_p_r_1_x=mesh_p_sec_2_x+length(r_1_x)-1;
mesh_p_r_1_y=mesh_p_sec_1_y+length(r_1_y);

%region 2

r_2_x=(x_s_2:dx:x_e_2)';
r_2_y=(y_s_2:dy:y_e_2)';

mesh_p_r_2_x=mesh_p_r_1_x+length(r_2_x)-1;
mesh_p_r_2_y=mesh_p_sec_1_y+length(r_2_y)-1;

%region 3

r_3_x=(x_s_3:dx:x_e_3)';
r_3_y=(y_s_3:dy:y_e_3)';

mesh_p_r_3_x=mesh_p_r_1_x+length(r_3_x)-1;
mesh_p_r_3_y=mesh_p_r_1_y+length(r_3_y)-1;

%region 4

r_4_x=(x_s_4:dx:x_e_4)';
r_4_y=(y_s_4:dy:y_e_4)';

mesh_p_r_4_x=mesh_p_sec_2_x+length(r_4_x)-1;
mesh_p_r_4_y=mesh_p_r_1_y+length(r_4_y)-1;

%region 5 section 3 and 4

% region 5 section 1 and 2
sec_3_x=(0:dx:X)';
sec_3_y=(Y-18:dy:Y)';

mesh_p_sec_3_x=length(sec_3_x)-1;
mesh_p_sec_3_y=mesh_p_sec_2_y+length(sec_3_y)-1;

sec_4_x=(x_e_3:dx:X)';
sec_4_y=(18:dy:Y-18)';

mesh_p_sec_4_x=mesh_p_r_3_x+length(sec_4_x)-1;
mesh_p_sec_4_y=length(sec_4_y)-1+mesh_p_sec_1_y;


%% material property in different meshes









