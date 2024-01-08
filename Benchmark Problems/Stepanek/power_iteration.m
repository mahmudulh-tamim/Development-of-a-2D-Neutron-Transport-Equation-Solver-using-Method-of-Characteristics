
function [k_new, flux_new]=power_iteration()

N_a=32;

[X,Y,dx,dy,sigma_t,sigma_s,nu_sigma_f]=geometry_data_structuring();


[exponential_portion,s_len,sum_s_len,alt_azim_theta,fin_d]=ray_tracing_polar_included(X,Y,dx,dy,N_a,sigma_t);




%given data
tol=10^(-7);





x=(0:dx:X)';
y=(0:dx:Y)';
n_x=length(x);
n_y=length(y);


%%

flux_old=(1/sum(sum(dx*dy*nu_sigma_f)))*ones(n_y-1,n_x-1);
k_old=1;

flux_new_half=source_iteration(flux_old, k_old,exponential_portion,s_len,sum_s_len,alt_azim_theta,fin_d,X,Y,dx,dy,N_a,sigma_t,sigma_s,nu_sigma_f);
k_new=k_old*sum(sum(nu_sigma_f*dx*dy.*flux_new_half))/sum(sum(nu_sigma_f*dx*dy.*flux_old));

flux_new=k_old/k_new*flux_new_half;

iteration=1;

while abs(max(max(flux_old-flux_new)))>tol %abs((k_new-k_old))>tol
    flux_old=flux_new;
    k_old=k_new;
    flux_new_half=source_iteration(flux_old, k_old,exponential_portion,s_len,sum_s_len,alt_azim_theta,fin_d,X,Y,dx,dy,N_a,sigma_t,sigma_s,nu_sigma_f);

    k_new=k_old*sum(sum(nu_sigma_f*dx*dy.*flux_new_half))/sum(sum(nu_sigma_f*dx*dy.*flux_old));

    flux_new=k_old/k_new*flux_new_half;

    iteration=1+iteration
end

mesh_centre_x=(dx/2:dx:X)';
mesh_centre_y=(dy/2:dy:Y)';
figure(2)

[mesh_x, mesh_y]=meshgrid(mesh_centre_x,mesh_centre_y);

surf(mesh_x,mesh_y,flux_new);
xlabel("X ordinate");
ylabel("Y ordinate");
zlabel("Flux");
set(gca,'zscale','log')
iteration