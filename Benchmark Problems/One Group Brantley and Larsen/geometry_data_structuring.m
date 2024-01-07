function [X,Y,dx,dy,sigma_t,sigma_s,nu_sigma_f]=geometry_data_structuring()

%given data

%% material data
%total cross section

sigma_t_moderator=1;
sigma_t_fuel=1.5;

%scattering cross section
sigma_s_moderator=0.93;
sigma_s_fuel=1.35;

%fission data
nu_sigma_f_moderator=0;
nu_sigma_f_fuel=0.24;

%geometry data
X=10;
Y=10;

%horizontal region

horizontal_moderator_start_x=0;
horizontal_moderator_start_y=0;

horizontal_moderator_end_x=10;
horizontal_moderator_end_y=1;

%region 1 moderator
vertical_moderator_region1_start_x=0;
vertical_moderator_region1_start_y=1;

vertical_moderator_region1_end_x=1;
vertical_moderator_region1_end_y=10;

%region 1 fuel

vertical_fuel_region1_start_x=vertical_moderator_region1_end_x;
vertical_fuel_region1_start_y=1;

vertical_fuel_region1_end_x=vertical_fuel_region1_start_x+1;
vertical_fuel_region1_end_y=10;

%region 2 moderator

vertical_moderator_region2_start_x=vertical_fuel_region1_end_x;
vertical_moderator_region2_start_y=1;

vertical_moderator_region2_end_x=vertical_moderator_region2_start_x+2;
vertical_moderator_region2_end_y=10;


%region 2 fuel

vertical_fuel_region2_start_x=vertical_moderator_region2_end_x;
vertical_fuel_region2_start_y=1;

vertical_fuel_region2_end_x=vertical_fuel_region2_start_x+1;
vertical_fuel_region2_end_y=10;

%region 3 moderator

vertical_moderator_region3_start_x=vertical_fuel_region2_end_x;
vertical_moderator_region3_start_y=1;

vertical_moderator_region3_end_x=vertical_moderator_region3_start_x+2;
vertical_moderator_region3_end_y=10;


%region 3 fuel

vertical_fuel_region3_start_x=vertical_moderator_region3_end_x;
vertical_fuel_region3_start_y=1;

vertical_fuel_region3_end_x=vertical_fuel_region3_start_x+1;
vertical_fuel_region3_end_y=10;

%region 4 moderator

vertical_moderator_region4_start_x=vertical_fuel_region3_end_x;
vertical_moderator_region4_start_y=1;

vertical_moderator_region4_end_x=vertical_moderator_region4_start_x+2;
vertical_moderator_region4_end_y=10;


%data structuring in matrix

dx=0.1; % integer multiple of dx should be equal to 1,2 and 10
dy=0.1; % integer multiple of dy should be equal to 1,9 and 10

x=(0:dx:X)';
y=(0:dx:Y)';

mesh_count_x=length(x)-1;
mesh_count_y=length(y)-1;

%initializtion of mesh properties
sigma_t=zeros(mesh_count_y,mesh_count_x);
sigma_s=zeros(mesh_count_y,mesh_count_x);
nu_sigma_f=zeros(mesh_count_y,mesh_count_x);

%% mesh count in different regions and material assignment

%horizontal region

hor_x=(horizontal_moderator_start_x:dx:horizontal_moderator_end_x)';
hor_y=(horizontal_moderator_start_y:dy:horizontal_moderator_end_y)';
hor_mesh_count_x=length(hor_x)-1;
hor_mesh_count_y=length(hor_y)-1;

sigma_t(1:hor_mesh_count_y,1:hor_mesh_count_x)=sigma_t_moderator;
sigma_s(1:hor_mesh_count_y,1:hor_mesh_count_x)=sigma_s_moderator;
nu_sigma_f(1:hor_mesh_count_y,1:hor_mesh_count_x)=nu_sigma_f_moderator;

%vertical moderator region 1
verM_1_x=(vertical_moderator_region1_start_x:dx:vertical_moderator_region1_end_x)';
verM_1_y=(vertical_moderator_region1_start_y:dy:vertical_moderator_region1_end_y)';
verM_1_mesh_c_x=length(verM_1_x)-1;
verM_1_mesh_c_y=length(verM_1_y)-1;

sigma_t(hor_mesh_count_y+1:hor_mesh_count_y+verM_1_mesh_c_y,1:verM_1_mesh_c_x)=sigma_t_moderator;
sigma_s(hor_mesh_count_y+1:hor_mesh_count_y+verM_1_mesh_c_y,1:verM_1_mesh_c_x)=sigma_s_moderator;
nu_sigma_f(hor_mesh_count_y+1:hor_mesh_count_y+verM_1_mesh_c_y,1:verM_1_mesh_c_x)=nu_sigma_f_moderator;

%vertical fuel region 1

verF_1_x=(vertical_fuel_region1_start_x:dx:vertical_fuel_region1_end_x);
verF_1_y=(vertical_fuel_region1_start_y:dy:vertical_fuel_region1_end_y);

verF_1_mesh_c_x=length(verF_1_x)-1;
verF_1_mesh_c_y=length(verF_1_y)-1;

sigma_t(hor_mesh_count_y+1:hor_mesh_count_y+verF_1_mesh_c_y,verM_1_mesh_c_x+1:verM_1_mesh_c_x+verF_1_mesh_c_x)=sigma_t_fuel;
sigma_s(hor_mesh_count_y+1:hor_mesh_count_y+verF_1_mesh_c_y,verM_1_mesh_c_x+1:verM_1_mesh_c_x+verF_1_mesh_c_x)=sigma_s_fuel;
nu_sigma_f(hor_mesh_count_y+1:hor_mesh_count_y+verF_1_mesh_c_y,verM_1_mesh_c_x+1:verM_1_mesh_c_x+verF_1_mesh_c_x)=nu_sigma_f_fuel;

%vertical moderator region 2


verM_2_x=(vertical_moderator_region2_start_x:dx:vertical_moderator_region2_end_x)';
verM_2_y=(vertical_moderator_region2_start_y:dy:vertical_moderator_region2_end_y)';
verM_2_mesh_c_x=length(verM_2_x)-1;
verM_2_mesh_c_y=length(verM_2_y)-1;

sigma_t(hor_mesh_count_y+1:hor_mesh_count_y+verM_2_mesh_c_y,verM_1_mesh_c_x+verF_1_mesh_c_x+1:verM_1_mesh_c_x+verF_1_mesh_c_x+verM_2_mesh_c_x)=sigma_t_moderator;
sigma_s(hor_mesh_count_y+1:hor_mesh_count_y+verM_2_mesh_c_y,verM_1_mesh_c_x+verF_1_mesh_c_x+1:verM_1_mesh_c_x+verF_1_mesh_c_x+verM_2_mesh_c_x)=sigma_s_moderator;
nu_sigma_f(hor_mesh_count_y+1:hor_mesh_count_y+verM_2_mesh_c_y,verM_1_mesh_c_x+verF_1_mesh_c_x+1:verM_1_mesh_c_x+verF_1_mesh_c_x+verM_2_mesh_c_x)=nu_sigma_f_moderator;


%vertical fuel region 2

verF_2_x=(vertical_fuel_region2_start_x:dx:vertical_fuel_region2_end_x);
verF_2_y=(vertical_fuel_region2_start_y:dy:vertical_fuel_region2_end_y);

verF_2_mesh_c_x=length(verF_2_x)-1;
verF_2_mesh_c_y=length(verF_2_y)-1;

sigma_t(hor_mesh_count_y+1:hor_mesh_count_y+verF_2_mesh_c_y,verM_1_mesh_c_x+verF_1_mesh_c_x+verM_2_mesh_c_x+1:verM_1_mesh_c_x+verF_1_mesh_c_x+verM_2_mesh_c_x+verF_2_mesh_c_x)=sigma_t_fuel;
sigma_s(hor_mesh_count_y+1:hor_mesh_count_y+verF_2_mesh_c_y,verM_1_mesh_c_x+verF_1_mesh_c_x+verM_2_mesh_c_x+1:verM_1_mesh_c_x+verF_1_mesh_c_x+verM_2_mesh_c_x+verF_2_mesh_c_x)=sigma_s_fuel;
nu_sigma_f(hor_mesh_count_y+1:hor_mesh_count_y+verF_2_mesh_c_y,verM_1_mesh_c_x+verF_1_mesh_c_x+verM_2_mesh_c_x+1:verM_1_mesh_c_x+verF_1_mesh_c_x+verM_2_mesh_c_x+verF_2_mesh_c_x)=nu_sigma_f_fuel;




%vertical moderator region 3


verM_3_x=(vertical_moderator_region3_start_x:dx:vertical_moderator_region3_end_x)';
verM_3_y=(vertical_moderator_region3_start_y:dy:vertical_moderator_region3_end_y)';
verM_3_mesh_c_x=length(verM_3_x)-1;
verM_3_mesh_c_y=length(verM_3_y)-1;

sigma_t(hor_mesh_count_y+1:hor_mesh_count_y+verM_3_mesh_c_y,verM_1_mesh_c_x+verF_1_mesh_c_x+verM_2_mesh_c_x+verF_2_mesh_c_x+1:verM_1_mesh_c_x+verF_1_mesh_c_x+verM_2_mesh_c_x+verF_2_mesh_c_x+verM_3_mesh_c_x)=sigma_t_moderator;
sigma_s(hor_mesh_count_y+1:hor_mesh_count_y+verM_3_mesh_c_y,verM_1_mesh_c_x+verF_1_mesh_c_x+verM_2_mesh_c_x+verF_2_mesh_c_x+1:verM_1_mesh_c_x+verF_1_mesh_c_x+verM_2_mesh_c_x+verF_2_mesh_c_x+verM_3_mesh_c_x)=sigma_s_moderator;
nu_sigma_f(hor_mesh_count_y+1:hor_mesh_count_y+verM_3_mesh_c_y,verM_1_mesh_c_x+verF_1_mesh_c_x+verM_2_mesh_c_x+verF_2_mesh_c_x+1:verM_1_mesh_c_x+verF_1_mesh_c_x+verM_2_mesh_c_x+verF_2_mesh_c_x+verM_3_mesh_c_x)=nu_sigma_f_moderator;


%vertical fuel region 3

verF_3_x=(vertical_fuel_region3_start_x:dx:vertical_fuel_region3_end_x);
verF_3_y=(vertical_fuel_region3_start_y:dy:vertical_fuel_region3_end_y);

verF_3_mesh_c_x=length(verF_3_x)-1;
verF_3_mesh_c_y=length(verF_3_y)-1;


sigma_t(hor_mesh_count_y+1:hor_mesh_count_y+verF_3_mesh_c_y,verM_1_mesh_c_x+verF_1_mesh_c_x+verM_2_mesh_c_x+verF_2_mesh_c_x+verM_3_mesh_c_x+1:verM_1_mesh_c_x+verF_1_mesh_c_x+verM_2_mesh_c_x+verF_2_mesh_c_x+verM_3_mesh_c_x+verF_3_mesh_c_x)=sigma_t_fuel;
sigma_s(hor_mesh_count_y+1:hor_mesh_count_y+verF_3_mesh_c_y,verM_1_mesh_c_x+verF_1_mesh_c_x+verM_2_mesh_c_x+verF_2_mesh_c_x+verM_3_mesh_c_x+1:verM_1_mesh_c_x+verF_1_mesh_c_x+verM_2_mesh_c_x+verF_2_mesh_c_x+verM_3_mesh_c_x+verF_3_mesh_c_x)=sigma_s_fuel;
nu_sigma_f(hor_mesh_count_y+1:hor_mesh_count_y+verF_3_mesh_c_y,verM_1_mesh_c_x+verF_1_mesh_c_x+verM_2_mesh_c_x+verF_2_mesh_c_x+verM_3_mesh_c_x+1:verM_1_mesh_c_x+verF_1_mesh_c_x+verM_2_mesh_c_x+verF_2_mesh_c_x+verM_3_mesh_c_x+verF_3_mesh_c_x)=nu_sigma_f_fuel;

%vertical moderator region 4


verM_4_x=(vertical_moderator_region4_start_x:dx:vertical_moderator_region4_end_x)';
verM_4_y=(vertical_moderator_region4_start_y:dy:vertical_moderator_region4_end_y)';
verM_4_mesh_c_x=length(verM_4_x)-1;
verM_4_mesh_c_y=length(verM_4_y)-1;


sigma_t(hor_mesh_count_y+1:hor_mesh_count_y+verM_4_mesh_c_y,verM_1_mesh_c_x+verF_1_mesh_c_x+verM_2_mesh_c_x+verF_2_mesh_c_x+verM_3_mesh_c_x+verF_3_mesh_c_x+1:verM_1_mesh_c_x+verF_1_mesh_c_x+verM_2_mesh_c_x+verF_2_mesh_c_x+verM_3_mesh_c_x+verF_3_mesh_c_x+verM_4_mesh_c_x)=sigma_t_moderator;
sigma_s(hor_mesh_count_y+1:hor_mesh_count_y+verM_4_mesh_c_y,verM_1_mesh_c_x+verF_1_mesh_c_x+verM_2_mesh_c_x+verF_2_mesh_c_x+verM_3_mesh_c_x+verF_3_mesh_c_x+1:verM_1_mesh_c_x+verF_1_mesh_c_x+verM_2_mesh_c_x+verF_2_mesh_c_x+verM_3_mesh_c_x+verF_3_mesh_c_x+verM_4_mesh_c_x)=sigma_s_moderator;
nu_sigma_f(hor_mesh_count_y+1:hor_mesh_count_y+verM_4_mesh_c_y,verM_1_mesh_c_x+verF_1_mesh_c_x+verM_2_mesh_c_x+verF_2_mesh_c_x+verM_3_mesh_c_x+verF_3_mesh_c_x+1:verM_1_mesh_c_x+verF_1_mesh_c_x+verM_2_mesh_c_x+verF_2_mesh_c_x+verM_3_mesh_c_x+verF_3_mesh_c_x+verM_4_mesh_c_x)=nu_sigma_f_moderator;



