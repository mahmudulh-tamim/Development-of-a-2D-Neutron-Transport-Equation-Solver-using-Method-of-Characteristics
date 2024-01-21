
%function [k_new, flux_new]=power_iteration()
[X,Y,dx,dy,sigma_t,sigma_s,nu_sigma_f]=geometry_data_structuring();

N_a=32;

[F_1,F_2,G_1,G_2,H,tau,ksi,x_c_t,y_c_t,X_i_c,Y_i_c,s_len,sum_s_len,adj_len,alt_azim_theta,fin_d, mesh_center_abscissa_number,mesh_center_ordinate_number, total_rays,M,c_i_xx, c_i_yy,c_i_xy]=ray_tracing_for_other_data(X,Y,dx,dy,N_a,sigma_t);

uuf=max(max(total_rays));
polar_discretization_number=3;
psi_bound=zeros(N_a,polar_discretization_number,uuf);


%given data
tol=10^(-6);



%spatial discretization



x=(0:dx:X)';
y=(0:dx:Y)';
n_x=length(x);
n_y=length(y);


%%

flux_old=(1/sum(sum(ones(n_x-1,n_y-1)*dx*dy.*nu_sigma_f)))*ones(n_x-1,n_y-1);
x_moment_flux=ones(mesh_center_abscissa_number,mesh_center_ordinate_number);
y_moment_flux=ones(mesh_center_abscissa_number,mesh_center_ordinate_number);
k_old=1;

[flux_new_half,x_moment_flux,y_moment_flux,psi_bound]=source_iteration(flux_old,x_moment_flux,y_moment_flux,k_old,F_1,F_2,G_1,G_2,H,tau,ksi,x_c_t,y_c_t,X_i_c,Y_i_c,s_len,sum_s_len,adj_len,alt_azim_theta,fin_d, mesh_center_abscissa_number,mesh_center_ordinate_number, total_rays,M,c_i_xx, c_i_yy,c_i_xy,N_a,X,Y,dx,dy,sigma_t,sigma_s,nu_sigma_f,psi_bound);
k_new=k_old*sum(sum(nu_sigma_f*dx*dy.*flux_new_half))/sum(sum(nu_sigma_f.*dx*dy.*flux_old));

flux_new=k_old/k_new*flux_new_half;

iteration=1;

while abs(max(max(flux_old-flux_new)))>tol %abs((k_new-k_old))>tol
    flux_old=flux_new;
    k_old=k_new;
    [flux_new_half,x_moment_flux,y_moment_flux,psi_bound]=source_iteration(flux_old,x_moment_flux,y_moment_flux,k_old,F_1,F_2,G_1,G_2,H,tau,ksi,x_c_t,y_c_t,X_i_c,Y_i_c,s_len,sum_s_len,adj_len,alt_azim_theta,fin_d, mesh_center_abscissa_number,mesh_center_ordinate_number, total_rays,M,c_i_xx, c_i_yy,c_i_xy,N_a,X,Y,dx,dy,sigma_t,sigma_s,nu_sigma_f,psi_bound);
    k_new=k_old*sum(sum(nu_sigma_f.*dx*dy.*flux_new_half))/sum(sum(nu_sigma_f.*dx*dy.*flux_old));

    flux_new=k_old/k_new*flux_new_half;

    

    iteration=1+iteration
end

mesh_centre_x=(dx/2:dx:X)';
mesh_centre_y=(dx/2:dy:Y)';
figure(50)
surfc(mesh_centre_x, mesh_centre_y,  flux_new');
xlabel("X ordinate");
ylabel("Y ordinate");
zlabel("Flux");

colorbar
iteration