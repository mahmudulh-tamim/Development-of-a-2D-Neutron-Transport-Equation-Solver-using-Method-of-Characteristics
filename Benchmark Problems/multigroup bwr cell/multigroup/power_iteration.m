
%function [k_new, flux_new]=power_iteration()

N_a=40;
tot_g=2;

[X,Y,dx,dy,sigma_t,sigma_s,nu_sigma_f,chi]=geometry_data_structuring();


[exponential_portion1,s_len,sum_s_len,alt_azim_theta,fin_d, mesh_center_ordinate_number,mesh_center_abscissa_number,total_rays]=ray_tracing_polar_included(X,Y,dx,dy,N_a,sigma_t(:,:,1));
[exponential_portion2,s_len,sum_s_len,alt_azim_theta,fin_d, mesh_center_ordinate_number,mesh_center_abscissa_number,total_rays]=ray_tracing_polar_included(X,Y,dx,dy,N_a,sigma_t(:,:,2));

jkaj=size(exponential_portion1);
size_exponential_portion=vertcat(jkaj',2);
exponential_portion=zeros(size_exponential_portion');
exponential_portion(:,:,:,:,:,1)=exponential_portion1;
exponential_portion(:,:,:,:,:,2)=exponential_portion2;


mesh_center_x=(dx/2:dx:X)';
mesh_center_y=(dy/2:dy:Y)';


%given data
tol=10^(-10);








%%

flux_old=ones(mesh_center_ordinate_number,mesh_center_abscissa_number,2);
k_old=1;
old_q_f=0;
for i=1:tot_g
    old_q_f=old_q_f+sum(sum(dx*dy*nu_sigma_f(:,:,i).*flux_old(:,:,i)));
end

flux_new=group_flux(flux_old, k_old,exponential_portion,s_len,sum_s_len,alt_azim_theta,fin_d,X,Y,dx,dy,N_a,sigma_t,sigma_s,nu_sigma_f,chi, mesh_center_ordinate_number,mesh_center_abscissa_number, total_rays);
new_q_f=0;
for i=1:tot_g
    new_q_f=new_q_f+sum(sum(dx*dy*nu_sigma_f(:,:,i).*flux_new(:,:,i)));
end
k_new=k_old*new_q_f/old_q_f;



iteration_power=1;

while max(max(max(abs(flux_new-flux_old))))>tol
    flux_old=flux_new;
    k_old=k_new;
    old_q_f=new_q_f;
   
    flux_new=group_flux(flux_old, k_old,exponential_portion,s_len,sum_s_len,alt_azim_theta,fin_d,X,Y,dx,dy,N_a,sigma_t,sigma_s,nu_sigma_f,chi, mesh_center_ordinate_number,mesh_center_abscissa_number, total_rays);
    new_q_f=0;
    for i=1:tot_g
        new_q_f=new_q_f+sum(sum(dx*dy*nu_sigma_f(:,:,i).*flux_new(:,:,i)));
    end

    k_new=k_old*new_q_f/old_q_f;

    iteration_power=iteration_power+1;
end

figure(1)
%[mesh_x, mesh_y]=meshgrid(mesh_center_abscissa_number,mesh_center_ordinate_number);
mesh_x=mesh_center_x;
mesh_y=mesh_center_y;
surfc(mesh_y,mesh_x,flux_new(:,:,1));

xlabel("X ordinate");
ylabel("Y ordinate");
zlabel("Flux");
title("Fast Flux");
colorbar
figure(2)
%[mesh_x, mesh_y]=meshgrid(mesh_center_abscissa_number,mesh_center_ordinate_number);
surfc(mesh_y,mesh_x,flux_new(:,:,2));

xlabel("X ordinate");
ylabel("Y ordinate");  
zlabel("Flux");
title("Thermal Flux");
colorbar
iteration_power