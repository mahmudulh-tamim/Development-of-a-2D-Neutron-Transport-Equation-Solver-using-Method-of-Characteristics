function flux_new=group_flux(flux_old, k_old,exponential_portion,s_len,sum_s_len,alt_azim_theta,fin_d,X,Y,dx,dy,N_a,sigma_t,sigma_s,nu_sigma_f,chi, mesh_center_ordinate_number,mesh_center_abscissa_number, total_rays)


tot_g=2;
flux_new=zeros(mesh_center_ordinate_number,mesh_center_abscissa_number,tot_g);

for g=tot_g:-1:1
    flux_new(:,:,g)=source_iteration(flux_old, k_old,exponential_portion,s_len,sum_s_len,alt_azim_theta,fin_d,X,Y,dx,dy,N_a,sigma_t,sigma_s,nu_sigma_f,chi,mesh_center_ordinate_number,mesh_center_abscissa_number,total_rays,g);
    flux_old(:,:,g)=flux_new(:,:,g);
end
