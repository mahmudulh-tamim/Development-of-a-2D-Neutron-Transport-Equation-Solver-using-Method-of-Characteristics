
%function [k_new, flux_new]=power_iteration()
X=4;
Y=4;
dx=0.5;
dy=0.5;
N_a=32;
sigma_t=1;
sigma_s=0.7;
nu_sigma_f=0.39;

[F_1,F_2,G_1,G_2,H,tau,ksi,x_c_t,y_c_t,X_i_c,Y_i_c,s_len,sum_s_len,adj_len,alt_azim_theta,fin_d, mesh_center_abscissa_number,mesh_center_ordinate_number, total_rays,M,c_i_xx, c_i_yy,c_i_xy]=ray_tracing_for_other_data(X,Y,dx,dy,N_a,sigma_t);



%given data
tol=10^(-6);



%spatial discretization



x=(0:dx:X)';
y=(0:dx:Y)';
n_x=length(x);
n_y=length(y);


%%

flux_old=(1/sum(sum(ones(n_x-1,n_y-1)*dx*dy*nu_sigma_f)))*ones(n_x-1,n_y-1);
x_moment_flux=ones(mesh_center_abscissa_number,mesh_center_ordinate_number);
y_moment_flux=ones(mesh_center_abscissa_number,mesh_center_ordinate_number);
k_old=1;

[flux_new_half,x_moment_flux,y_moment_flux]=source_iteration(flux_old,x_moment_flux,y_moment_flux,k_old,F_1,F_2,G_1,G_2,H,tau,ksi,x_c_t,y_c_t,X_i_c,Y_i_c,s_len,sum_s_len,adj_len,alt_azim_theta,fin_d, mesh_center_abscissa_number,mesh_center_ordinate_number, total_rays,M,c_i_xx, c_i_yy,c_i_xy,N_a);
k_new=k_old*sum(sum(nu_sigma_f*dx*dy*flux_new_half))/sum(sum(nu_sigma_f*dx*dy*flux_old));

flux_new=k_old/k_new*flux_new_half;

iteration=1;

while abs(max(max(flux_old-flux_new)))>tol %abs((k_new-k_old))>tol
    flux_old=flux_new;
    k_old=k_new;
    [flux_new_half,x_moment_flux,y_moment_flux]=source_iteration(flux_old,x_moment_flux,y_moment_flux,k_old,F_1,F_2,G_1,G_2,H,tau,ksi,x_c_t,y_c_t,X_i_c,Y_i_c,s_len,sum_s_len,adj_len,alt_azim_theta,fin_d, mesh_center_abscissa_number,mesh_center_ordinate_number, total_rays,M,c_i_xx, c_i_yy,c_i_xy,N_a);
    k_new=k_old*sum(sum(nu_sigma_f*dx*dy*flux_new_half))/sum(sum(nu_sigma_f*dx*dy*flux_old));

    flux_new=k_old/k_new*flux_new_half;

    

    iteration=1+iteration
end

mesh_centre_x=(dx/2:dx:X)';
mesh_centre_y=(dx/2:dy:Y)';
figure(50)
mesh(mesh_centre_x, mesh_centre_y,  flux_new);
xlabel("X ordinate");
ylabel("Y ordinate");
zlabel("Flux");

iteration