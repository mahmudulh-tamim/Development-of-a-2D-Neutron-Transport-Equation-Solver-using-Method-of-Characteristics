function flux_new=source_iteration(flux_old, k_old,exponential_portion,s_len,sum_s_len,altered_azimuthal_direction_theta,fin_d)

%given data
tol=10^(-5);


sigma_t=1;
sigma_s=0.7;
nu_sigma_f=0.39;

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
fission_term= (1/k_old)*nu_sigma_f*flux_old;

scattering_term=sigma_s*flux_old;

source_term=(fission_term+scattering_term)/(4*pi);

flux_new=transport_sweep(source_term,exponential_portion,s_len,sum_s_len,altered_azimuthal_direction_theta,fin_d);

iteration=1;

while max(max(abs(flux_new-flux_old)))>tol
    flux_old=flux_new;
    scattering_term=sigma_s*flux_old;

    source_term=(fission_term+scattering_term)/(4*pi);

    flux_new=transport_sweep(source_term,exponential_portion,s_len,sum_s_len,altered_azimuthal_direction_theta,fin_d);

    iteration=1+iteration;
end
