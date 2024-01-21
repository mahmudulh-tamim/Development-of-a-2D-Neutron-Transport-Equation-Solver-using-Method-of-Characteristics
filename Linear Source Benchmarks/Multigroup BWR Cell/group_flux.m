function [flux_new,new_x_moment_flux,new_y_moment_flux,new_psi_bound]=group_flux(flux_old,x_moment_flux,y_moment_flux,k_old,F_1,F_2,G_1,G_2,H,tau,ksi,x_c_t,y_c_t,X_i_c,Y_i_c,s_len,sum_s_len,adj_len,alt_azim_theta,fin_d, mesh_center_abscissa_number,mesh_center_ordinate_number, total_rays,M,c_i_xx, c_i_yy,c_i_xy,N_a,X,Y,dx,dy,sigma_t,sigma_s,nu_sigma_f,chi,psi_bound)


tot_g=2;
flux_new=zeros(mesh_center_abscissa_number,mesh_center_ordinate_number,tot_g);
new_x_moment_flux=zeros(mesh_center_abscissa_number,mesh_center_ordinate_number,tot_g);
new_y_moment_flux=zeros(mesh_center_abscissa_number,mesh_center_ordinate_number,tot_g);

uuf=max(max(total_rays));
polar_discretization_number=3;
new_psi_bound=zeros(N_a,polar_discretization_number,uuf,tot_g);

for g=tot_g:-1:1
    [flux_new(:,:,g),new_x_moment_flux(:,:,g),new_y_moment_flux(:,:,g),new_psi_bound(:,:,:,g)]=source_iteration(flux_old,x_moment_flux,y_moment_flux,k_old,F_1(:,:,:,:,:,g),F_2(:,:,:,:,:,g),G_1(:,:,:,:,:,g),G_2(:,:,:,:,:,g),H(:,:,:,:,:,g),tau(:,:,:,:,:,g),ksi,x_c_t,y_c_t,X_i_c,Y_i_c,s_len,sum_s_len,adj_len,alt_azim_theta,fin_d, mesh_center_abscissa_number,mesh_center_ordinate_number, total_rays,M,c_i_xx(:,:,g), c_i_yy(:,:,g),c_i_xy(:,:,g),N_a,X,Y,dx,dy,sigma_t(:,:,g),sigma_s,nu_sigma_f,chi,psi_bound(:,:,:,g),g);
    
end
