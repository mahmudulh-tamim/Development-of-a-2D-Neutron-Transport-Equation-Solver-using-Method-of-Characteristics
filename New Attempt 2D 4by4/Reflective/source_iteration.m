function flux_new=source_iteration(flux_old, k_old,exponential_portion,s_len,sum_s_len,alt_azim_theta,fin_d,N_a,total_rays)

%given data
tol=10^(-10);


sigma_t=1;
sigma_s=0.7;
nu_sigma_f=0.39;

%spatial discretization

X=4;
Y=4;

dx=0.1;
dy=0.1;

mesh_center_x=(dx/2:dx:X)';
mesh_center_y=(dy/2:dy:Y)';
mesh_center_abscissa_number=length(mesh_center_x);
mesh_center_ordinate_number=length(mesh_center_y);



%%
fission_term= (1/k_old)*nu_sigma_f*flux_old;

scattering_term=sigma_s*flux_old;

source_term=(fission_term+scattering_term)/(4*pi);

flux_new=reflective(source_term,exponential_portion,s_len,sum_s_len,alt_azim_theta,fin_d,X,Y,dx,dy,N_a,sigma_t,total_rays);

iteration=1;
%{
while max(max(abs(flux_new-flux_old)))>tol
    flux_old=flux_new;
    scattering_term=sigma_s*flux_old;

    source_term=(fission_term+scattering_term)/(4*pi);

    flux_new=reflective(source_term,exponential_portion,s_len,sum_s_len,alt_azim_theta,fin_d,X,Y,dx,dy,N_a,sigma_t,total_rays);

    iteration=1+iteration;
end
%}
