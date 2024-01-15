function flux_gth_group=source_iteration(flux_old, k_old,exponential_portion,s_len,sum_s_len,alt_azim_theta,fin_d,X,Y,dx,dy,N_a,sigma_t,sigma_s,nu_sigma_f,chi,mesh_center_ordinate_number,mesh_center_abscissa_number,total_rays,g)

%given data
tol=10^(-10);
tot_g=2;

mesh_center_x=(dx/2:dx:X)';
mesh_center_y=(dy/2:dy:Y)';


sig_t=sigma_t(:,:,g);

flux_new=flux_old;

%%
fission_term=zeros(mesh_center_ordinate_number,mesh_center_abscissa_number);
for i=1:tot_g
    fission_term=fission_term+(1/k_old)*chi(:,:,g).*nu_sigma_f(:,:,i).*flux_old(:,:,i);
end
scattering_term_ind=zeros(mesh_center_ordinate_number,mesh_center_abscissa_number);
for i=1:tot_g
    if(i~=g) 
        scattering_term_ind=scattering_term_ind+sigma_s(:,:,i,g).*flux_old(:,:,i);
    end
end
scattering_term_self=sigma_s(:,:,g,g).*flux_old(:,:,g);

source_term=(fission_term+scattering_term_ind+scattering_term_self)/(4*pi);


flux_new(:,:,g)=reflective(source_term,exponential_portion,s_len,sum_s_len,alt_azim_theta,fin_d,X,Y,dx,dy,N_a,sig_t,total_rays);



iteration_source=1;

while max(max(abs(flux_new(:,:,g)-flux_old(:,:,g))))>tol
    flux_old(:,:,g)=flux_new(:,:,g);
    
    scattering_term_self=sigma_s(:,:,g,g).*flux_old(:,:,g);
    source_term=(fission_term+scattering_term_ind+scattering_term_self)/(4*pi);

    flux_new(:,:,g)=reflective(source_term,exponential_portion,s_len,sum_s_len,alt_azim_theta,fin_d,X,Y,dx,dy,N_a,sig_t,total_rays);

    
    iteration_source=1+iteration_source;
end

flux_gth_group=flux_new(:,:,g);
iteration_source
