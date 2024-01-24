function [flux_new, x_moment_flux, y_moment_flux,psi_bound]=source_iteration(flux_old,x_moment_flux,y_moment_flux,k_old,F_1,F_2,H,ksi,x_c_t,y_c_t,X_i_c,Y_i_c,s_len,adj_len,alt_azim_theta,fin_d, mesh_center_abscissa_number,mesh_center_ordinate_number, total_rays,M,c_i_xx, c_i_yy,c_i_xy,N_a,X,Y,dx,dy,sigma_t,sigma_s,nu_sigma_f,psi_bound)
%given data
tol=10^(-10);

%spatial discretization

mesh_center_x=(dx/2:dx:X)';
mesh_center_y=(dy/2:dy:Y)';
mesh_center_abscissa_number=length(mesh_center_x);
mesh_center_ordinate_number=length(mesh_center_y);

%% angular discretization

%polar discretization

mu=[0.932954;0.537707;0.166648];
w=2*[0.670148;0.283619;0.046233];
polar_discretization_number=size(mu,1);

%azimuthal discretization

del_theta=2*pi/N_a;
theta=(0:del_theta:2*pi)';
azimuthal_direction_theta= 0.5*(theta(1:end-1,1)+theta(2:end,1));
azimuthal_discretization_number=size(azimuthal_direction_theta,1);


q_bar_i_fission=nu_sigma_f.*flux_old/(k_old);
q_bar_i=(sigma_s.*flux_old+q_bar_i_fission);
Q_i_x_fission=nu_sigma_f.*x_moment_flux/k_old;
Q_i_x=(sigma_s.*x_moment_flux+Q_i_x_fission);
Q_i_y_fission=nu_sigma_f.*y_moment_flux/k_old;
Q_i_y=(sigma_s.*y_moment_flux+Q_i_y_fission);

%q_bar_i_f=zeros(mesh_center_abscissa_number,mesh_center_ordinate_number);
q_i_x=zeros(mesh_center_abscissa_number,mesh_center_ordinate_number);
q_i_y=zeros(mesh_center_abscissa_number,mesh_center_ordinate_number);

for j=1:mesh_center_ordinate_number
    for i=1:mesh_center_abscissa_number
        temp_M=M(:,:,i,j);
        solved=temp_M\[q_bar_i(i,j); Q_i_x(i,j); Q_i_y(i,j)];
        
        q_bar_i(i,j)=solved(1,1);
        q_i_x(i,j)=solved(2,1);
        q_i_y(i,j)=solved(3,1);
    end
end






%%


[flux_new,x_moment_flux,y_moment_flux,psi_bound]=transport_sweep(q_i_x,q_i_y,q_bar_i,c_i_xx,c_i_yy,c_i_xy,F_1,F_2,H,adj_len,s_len,ksi,alt_azim_theta,fin_d,x_c_t,y_c_t,X_i_c,Y_i_c,X,Y,dx,dy,N_a,sigma_t,psi_bound);

iteration=1;

%{
while max(max(abs(flux_new-flux_old)))>tol
    flux_old=flux_new;
    q_bar_i=sigma_s*flux_old+q_bar_i_fission;
    Q_i_x=sigma_s*x_moment_flux+Q_i_x_fission;
    Q_i_y=sigma_s*y_moment_flux+Q_i_y_fission;
    for j=1:mesh_center_ordinate_number
        for i=1:mesh_center_abscissa_number
            temp_M=M(:,:,i,j);
            solved=temp_M\[q_bar_i(i,j); Q_i_x(i,j); Q_i_y(i,j)];
            
            q_bar_i(i,j)=solved(1,1);
            q_i_x(i,j)=solved(2,1);
            q_i_y(i,j)=solved(3,1);
        end
    end

    [flux_new,x_moment_flux,y_moment_flux]=transport_sweep(q_i_x,q_i_y,q_bar_i,c_i_xx,c_i_yy,c_i_xy,F_1,F_2,H,adj_len,s_len,ksi,alt_azim_theta,fin_d,x_c_t,y_c_t,X_i_c,Y_i_c,X,Y,dx,dy,N_a,sigma_t,psi_bound);


    iteration=1+iteration;
end
%}