%function [X,Y,dx,dy,sigma_t,sigma_s,nu_sigma_f,chi]=geometry_data_structuring()
%group number
g=2;
%region fuel
%group 1
chi_f_1=1;
sigma_t_f_1=0.196647;
sigma_s_f_1=[0.178, 0.001089]';

nu_sigma_f_f_1=0.006203;

%group 2
chi_f_2=0;
sigma_t_f_2=0.596159;
sigma_s_f_2=[0.01002, 0.5255]';
nu_sigma_f_f_2=0.1101;

%region moderator
%group 1
chi_m_1=0;
sigma_t_m_1=0.222064;
sigma_s_m_1=[0.1995, 0.001558]';

nu_sigma_f_m_1=0;

%group 2
chi_m_2=0;
sigma_t_m_2=0.887874;
sigma_s_m_2=[0.02118, 0.8783]';

nu_sigma_f_m_2=0;



%% geometry

X=8.9;
Y=8.9;

dx=0.05; %dx and dy should divide 18,25, 96,86 in integer
dy=0.05;

x=(0:dx:X)';
y=(0:dy:Y)';

mesh_count_x=length(x)-1;
mesh_count_y=length(y)-1;

%% initialization
chi=zeros(mesh_count_y,mesh_count_x,g);
sigma_t=zeros(mesh_count_y,mesh_count_x,g);
sigma_s=zeros(mesh_count_y,mesh_count_x, g,g);
nu_sigma_f=zeros(mesh_count_y,mesh_count_x,g);

%%

% region fuel
x_s_f=1.5;
y_s_f=1;

x_e_f=x_s_f+6.4;
y_e_f=y_s_f+6.4;

% region moderator

x_s_m=0;
y_s_m=0;

x_e_m=8.9;
y_e_m=8.9;

%% mesh geometry

% region moderator section 1 and 2
sec_1_x=(0:dx:X)';
sec_1_y=(0:dy:1)';

mesh_p_sec_1_x=length(sec_1_x)-1;
mesh_p_sec_1_y=length(sec_1_y)-1;

sec_2_x=(0:dx:1.5)';
sec_2_y=(1:dy:Y-1.5)';

mesh_p_sec_2_x=length(sec_2_x)-1;
mesh_p_sec_2_y=length(sec_2_y)-1+mesh_p_sec_1_y;

%region fuel

r_f_x=(x_s_f:dx:x_e_f);
r_f_y=(y_s_f:dy:y_e_f);

mesh_p_r_f_x=mesh_p_sec_2_x+length(r_f_x)-1;
mesh_p_r_f_y=mesh_p_sec_1_y+length(r_f_y)-1;



%region moderator section 3 and 4


sec_3_x=(0:dx:X)';
sec_3_y=(Y-1.5:dy:Y)';

mesh_p_sec_3_x=length(sec_3_x)-1;
mesh_p_sec_3_y=mesh_p_sec_2_y+length(sec_3_y)-1;

sec_4_x=(x_e_f:dx:X)';
sec_4_y=(1:dy:Y-1.5)';

mesh_p_sec_4_x=mesh_p_r_f_x+length(sec_4_x)-1;
mesh_p_sec_4_y=length(sec_4_y)-1+mesh_p_sec_1_y;


%% material property in different meshes

%region moderator


%section 1
%group 1
chi(1:mesh_p_sec_1_y,1:mesh_p_sec_1_x,1)=chi_m_1;
sigma_t(1:mesh_p_sec_1_y,1:mesh_p_sec_1_x,1)=sigma_t_m_1;
nu_sigma_f(1:mesh_p_sec_1_y,1:mesh_p_sec_1_x,1)=nu_sigma_f_m_1;
for i=1:g
    sigma_s(1:mesh_p_sec_1_y,1:mesh_p_sec_1_x,i,1)=sigma_s_m_1(i,1);
end

%group 2
chi(1:mesh_p_sec_1_y,1:mesh_p_sec_1_x,2)=chi_m_2;
sigma_t(1:mesh_p_sec_1_y,1:mesh_p_sec_1_x,2)=sigma_t_m_2;
nu_sigma_f(1:mesh_p_sec_1_y,1:mesh_p_sec_1_x,2)=nu_sigma_f_m_2;
for i=1:g
    sigma_s(1:mesh_p_sec_1_y,1:mesh_p_sec_1_x,i,2)=sigma_s_m_2(i,1);
end


%section 2
%group 1
chi(mesh_p_sec_1_y+1:mesh_p_sec_2_y,1:mesh_p_sec_2_x,1)=chi_m_1;
sigma_t(mesh_p_sec_1_y+1:mesh_p_sec_2_y,1:mesh_p_sec_2_x,1)=sigma_t_m_1;

nu_sigma_f(mesh_p_sec_1_y+1:mesh_p_sec_2_y,1:mesh_p_sec_2_x,1)=nu_sigma_f_m_1;
for i=1:g
    sigma_s(mesh_p_sec_1_y+1:mesh_p_sec_2_y,1:mesh_p_sec_2_x,i,1)=sigma_s_m_1(i,1);
end

%group 2
chi(mesh_p_sec_1_y+1:mesh_p_sec_2_y,1:mesh_p_sec_2_x,2)=chi_m_2;
sigma_t(mesh_p_sec_1_y+1:mesh_p_sec_2_y,1:mesh_p_sec_2_x,2)=sigma_t_m_2;

nu_sigma_f(mesh_p_sec_1_y+1:mesh_p_sec_2_y,1:mesh_p_sec_2_x,2)=nu_sigma_f_m_2;
for i=1:g
    sigma_s(mesh_p_sec_1_y+1:mesh_p_sec_2_y,1:mesh_p_sec_2_x,i,2)=sigma_s_m_2(i,1);
end



%section 3

%group 1
chi(mesh_p_sec_2_y+1:mesh_p_sec_3_y,1:mesh_p_sec_3_x,1)=chi_m_1;
sigma_t(mesh_p_sec_2_y+1:mesh_p_sec_3_y,1:mesh_p_sec_3_x,1)=sigma_t_m_1;

nu_sigma_f(mesh_p_sec_2_y+1:mesh_p_sec_3_y,1:mesh_p_sec_3_x,1)=nu_sigma_f_m_1;
for i=1:g
    sigma_s(mesh_p_sec_2_y+1:mesh_p_sec_3_y,1:mesh_p_sec_3_x,i,1)=sigma_s_m_1(i,1);
end


%group 2
chi(mesh_p_sec_2_y+1:mesh_p_sec_3_y,1:mesh_p_sec_3_x,2)=chi_m_2;
sigma_t(mesh_p_sec_2_y+1:mesh_p_sec_3_y,1:mesh_p_sec_3_x,2)=sigma_t_m_2;

nu_sigma_f(mesh_p_sec_2_y+1:mesh_p_sec_3_y,1:mesh_p_sec_3_x,2)=nu_sigma_f_m_2;
for i=1:g
    sigma_s(mesh_p_sec_2_y+1:mesh_p_sec_3_y,1:mesh_p_sec_3_x,i,2)=sigma_s_m_2(i,1);
end


%section 4
%group 1
chi(mesh_p_sec_1_y+1:mesh_p_sec_4_y,mesh_p_r_f_x+1:mesh_p_sec_4_x,1)=chi_m_1;
sigma_t(mesh_p_sec_1_y+1:mesh_p_sec_4_y,mesh_p_r_f_x+1:mesh_p_sec_4_x,1)=sigma_t_m_1;

nu_sigma_f(mesh_p_sec_1_y+1:mesh_p_sec_4_y,mesh_p_r_f_x+1:mesh_p_sec_4_x,1)=nu_sigma_f_m_1;
for i=1:g
    sigma_s(mesh_p_sec_1_y+1:mesh_p_sec_4_y,mesh_p_r_f_x+1:mesh_p_sec_4_x,i,1)=sigma_s_m_1(i,1);
end

%group 2
chi(mesh_p_sec_1_y+1:mesh_p_sec_4_y,mesh_p_r_f_x+1:mesh_p_sec_4_x,2)=chi_m_2;
sigma_t(mesh_p_sec_1_y+1:mesh_p_sec_4_y,mesh_p_r_f_x+1:mesh_p_sec_4_x,2)=sigma_t_m_2;

nu_sigma_f(mesh_p_sec_1_y+1:mesh_p_sec_4_y,mesh_p_r_f_x+1:mesh_p_sec_4_x,2)=nu_sigma_f_m_2;
for i=1:g
    sigma_s(mesh_p_sec_1_y+1:mesh_p_sec_4_y,mesh_p_r_f_x+1:mesh_p_sec_4_x,i,2)=sigma_s_m_2(i,1);
end



% region fuel
%group 1
chi(mesh_p_sec_1_y+1:mesh_p_r_f_y,mesh_p_sec_2_x+1:mesh_p_r_f_x,1)=chi_f_1;
sigma_t(mesh_p_sec_1_y+1:mesh_p_r_f_y,mesh_p_sec_2_x+1:mesh_p_r_f_x,1)=sigma_t_f_1;

nu_sigma_f(mesh_p_sec_1_y+1:mesh_p_r_f_y,mesh_p_sec_2_x+1:mesh_p_r_f_x,1)=nu_sigma_f_f_1;
for i=1:g
    sigma_s(mesh_p_sec_1_y+1:mesh_p_r_f_y,mesh_p_sec_2_x+1:mesh_p_r_f_x,i,1)=sigma_s_f_1(i,1);
end

%group 2
chi(mesh_p_sec_1_y+1:mesh_p_r_f_y,mesh_p_sec_2_x+1:mesh_p_r_f_x,2)=chi_f_2;
sigma_t(mesh_p_sec_1_y+1:mesh_p_r_f_y,mesh_p_sec_2_x+1:mesh_p_r_f_x,2)=sigma_t_f_2;

nu_sigma_f(mesh_p_sec_1_y+1:mesh_p_r_f_y,mesh_p_sec_2_x+1:mesh_p_r_f_x,2)=nu_sigma_f_f_2;
for i=1:g
    sigma_s(mesh_p_sec_1_y+1:mesh_p_r_f_y,mesh_p_sec_2_x+1:mesh_p_r_f_x,i,2)=sigma_s_f_2(i,1);
end

