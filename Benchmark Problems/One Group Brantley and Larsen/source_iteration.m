function flux_new=source_iteration(flux_old, k_old,exponential_portion,s_len,sum_s_len,alt_azim_theta,fin_d,X,Y,dx,dy,N_a,sigma_t,sigma_s,nu_sigma_f,total_rays)

%given data
tol=10^(-7);






%%
fission_term= (1/k_old)*nu_sigma_f.*flux_old;

scattering_term=sigma_s.*flux_old;

source_term=(fission_term+scattering_term)/(4*pi);

flux_new=reflective(source_term,exponential_portion,s_len,sum_s_len,alt_azim_theta,fin_d,X,Y,dx,dy,N_a,sigma_t,total_rays);



iteration=1;
%{
while max(max(abs(flux_new-flux_old)))>tol
    flux_old=flux_new;
    scattering_term=sigma_s*flux_old;

    source_term=(fission_term+scattering_term)/(4*pi);

    flux_new=transport_sweep(source_term,exponential_portion,s_len,sum_s_len,alt_azim_theta,fin_d,X,Y,dx,dy,N_a,sigma_t);


    iteration=1+iteration;
end
iteration
%}