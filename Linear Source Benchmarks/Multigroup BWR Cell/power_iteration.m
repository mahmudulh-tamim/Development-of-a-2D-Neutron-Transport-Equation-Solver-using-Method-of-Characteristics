
%function [k_new, flux_new]=power_iteration()
[X,Y,dx,dy,mesh_count_x,mesh_count_y,sigma_t,sigma_s,nu_sigma_f,chi]=geometry_data_structuring();

N_a=16;
tot_g=2;

c_i_xx=zeros(mesh_count_x,mesh_count_y,2);
c_i_yy=zeros(mesh_count_x,mesh_count_y,2);
c_i_xy=zeros(mesh_count_x,mesh_count_y,2);

[F_1_1,F_2_1,G_1_1,G_2_1,H_1,tau_1,ksi,x_c_t,y_c_t,X_i_c,Y_i_c,s_len,sum_s_len,adj_len,alt_azim_theta,fin_d, mesh_center_abscissa_number,mesh_center_ordinate_number, total_rays,M,c_i_xx(:,:,1), c_i_yy(:,:,1),c_i_xy(:,:,1)]=ray_tracing_for_other_data(X,Y,dx,dy,N_a,sigma_t(:,:,1));
[F_1_2,F_2_2,G_1_2,G_2_2,H_2,tau_2,ksi,x_c_t,y_c_t,X_i_c,Y_i_c,s_len,sum_s_len,adj_len,alt_azim_theta,fin_d, mesh_center_abscissa_number,mesh_center_ordinate_number, total_rays,M,c_i_xx(:,:,2), c_i_yy(:,:,2),c_i_xy(:,:,2)]=ray_tracing_for_other_data(X,Y,dx,dy,N_a,sigma_t(:,:,2));
size_F=size(F_1_1);
size_F=vertcat(size_F',2);
F_1=zeros(size_F');
F_2=zeros(size_F');
G_1=zeros(size_F');
G_2=zeros(size_F');
H=zeros(size_F');
tau=zeros(size_F');
F_1(:,:,:,:,:,1)=F_1_1;
F_1(:,:,:,:,:,2)=F_1_2;
F_2(:,:,:,:,:,1)=F_2_1;
F_2(:,:,:,:,:,2)=F_2_2;
G_1(:,:,:,:,:,1)=G_1_1;
G_1(:,:,:,:,:,2)=G_1_2;
G_2(:,:,:,:,:,1)=G_2_1;
G_2(:,:,:,:,:,2)=G_2_2;
H(:,:,:,:,:,1)=H_1;
H(:,:,:,:,:,2)=H_2;
tau(:,:,:,:,:,1)=tau_1;
tau(:,:,:,:,:,2)=tau_2;

uuf=max(max(total_rays));
polar_discretization_number=3;
psi_bound=zeros(N_a,polar_discretization_number,uuf,2);


%given data
tol=10^(-6);



%spatial discretization



x=(0:dx:X)';
y=(0:dx:Y)';
n_x=length(x);
n_y=length(y);


%%

flux_old=ones(mesh_center_abscissa_number,mesh_center_ordinate_number,2);
x_moment_flux=ones(mesh_center_abscissa_number,mesh_center_ordinate_number,2);
y_moment_flux=ones(mesh_center_abscissa_number,mesh_center_ordinate_number,2);
k_old=1;

[flux_new_half,x_moment_flux,y_moment_flux,psi_bound]=group_flux(flux_old,x_moment_flux,y_moment_flux,k_old,F_1,F_2,G_1,G_2,H,tau,ksi,x_c_t,y_c_t,X_i_c,Y_i_c,s_len,sum_s_len,adj_len,alt_azim_theta,fin_d, mesh_center_abscissa_number,mesh_center_ordinate_number, total_rays,M,c_i_xx, c_i_yy,c_i_xy,N_a,X,Y,dx,dy,sigma_t,sigma_s,nu_sigma_f,chi,psi_bound);
old_q_f=0;
new_q_f=0;
for i=1:tot_g
    old_q_f=old_q_f+sum(sum(dx*dy*nu_sigma_f(:,:,i).*flux_old(:,:,i)));
    new_q_f=new_q_f+sum(sum(dx*dy*nu_sigma_f(:,:,i).*flux_new_half(:,:,i)));
end
k_new=k_old*new_q_f/old_q_f;

flux_new=k_old/k_new*flux_new_half;

iteration=1

while abs((k_new-k_old))>tol
    flux_old=flux_new;
    k_old=k_new;
    old_q_f=new_q_f;
    [flux_new_half,x_moment_flux,y_moment_flux,psi_bound]=group_flux(flux_old,x_moment_flux,y_moment_flux,k_old,F_1,F_2,G_1,G_2,H,tau,ksi,x_c_t,y_c_t,X_i_c,Y_i_c,s_len,sum_s_len,adj_len,alt_azim_theta,fin_d, mesh_center_abscissa_number,mesh_center_ordinate_number, total_rays,M,c_i_xx, c_i_yy,c_i_xy,N_a,X,Y,dx,dy,sigma_t,sigma_s,nu_sigma_f,chi,psi_bound);
    for i=1:tot_g
        new_q_f=new_q_f+sum(sum(dx*dy*nu_sigma_f(:,:,i).*flux_new_half(:,:,i)));
    end
    k_new=k_old*new_q_f/old_q_f;

    flux_new=k_old/k_new*flux_new_half;

    iteration=1+iteration
end

mesh_centre_x=(dx/2:dx:X)';
mesh_centre_y=(dx/2:dy:Y)';
figure(1)
surfc(mesh_centre_x, mesh_centre_y,  flux_new(:,:,1)');
xlabel("X ordinate");
ylabel("Y ordinate");
zlabel("Flux");
title("Fast Flux");
colorbar

figure(2)
surfc(mesh_centre_x, mesh_centre_y,  flux_new(:,:,2)');
xlabel("X ordinate");
ylabel("Y ordinate");
zlabel("Flux");
title("Thermal Flux");
colorbar
