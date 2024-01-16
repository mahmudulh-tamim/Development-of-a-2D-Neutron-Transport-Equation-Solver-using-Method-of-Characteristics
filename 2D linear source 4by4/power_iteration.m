
%function [k_new, flux_new]=power_iteration()
X=4;
Y=4;
dx=0.1;
dy=0.1;
N_a=64;
sigma_t=1;
sigma_s=0.7;
nu_sigma_f=0.39;

[exponential_portion,s_len,sum_s_len,alt_azim_theta,fin_d, mesh_center_ordinate_number,mesh_center_abscissa_number,total_rays]=ray_tracing_polar_included(X,Y,dx,dy,N_a,sigma_t);



%given data
tol=10^(-6);



%spatial discretization

X=4;
Y=4;

dx=0.1;
dy=0.1;

x=(0:dx:X)';
y=(0:dx:Y)';
n_x=length(x);
n_y=length(y);


%%

flux_old=(1/sum(sum(ones(n_x-1,n_y-1)*dx*dy*nu_sigma_f)))*ones(n_x-1,n_y-1);
k_old=1;

flux_new_half=source_iteration(flux_old, k_old,exponential_portion,s_len,sum_s_len,alt_azim_theta,fin_d,N_a,total_rays);
k_new=k_old*sum(sum(nu_sigma_f*dx*dy*flux_new_half))/sum(sum(nu_sigma_f*dx*dy*flux_old));

flux_new=k_old/k_new*flux_new_half;

iteration=1;

while abs(max(max(flux_old-flux_new)))>tol %abs((k_new-k_old))>tol
    flux_old=flux_new;
    k_old=k_new;
    flux_new_half=source_iteration(flux_old, k_old,exponential_portion,s_len,sum_s_len,alt_azim_theta,fin_d,N_a,total_rays);
    k_new=k_old*sum(sum(nu_sigma_f*dx*dy*flux_new_half))/sum(sum(nu_sigma_f*dx*dy*flux_old));

    flux_new=k_old/k_new*flux_new_half;

    iteration=1+iteration;
end
mesh_centre_x=(dx/2:dx:X)';
mesh_centre_y=(dx/2:dx:Y)';
figure(50)
mesh(mesh_centre_y, mesh_centre_x, flux_new);
xlabel("X ordinate");
ylabel("Y ordinate");
zlabel("Flux");

iteration