function [gth_flux, gth_x_moment_flux, gth_y_moment_flux,psi_bound]=source_iteration(flux_old,x_moment_flux,y_moment_flux,k_old,F_1,F_2,G_1,G_2,H,tau,ksi,x_c_t,y_c_t,X_i_c,Y_i_c,s_len,sum_s_len,adj_len,alt_azim_theta,fin_d, mesh_center_abscissa_number,mesh_center_ordinate_number, total_rays,M,c_i_xx, c_i_yy,c_i_xy,N_a,X,Y,dx,dy,sigma_t,sigma_s,nu_sigma_f,chi,psi_bound,g)
%given data
tol=10^(-7);

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

q_bar_i_fission=zeros(mesh_center_abscissa_number,mesh_center_ordinate_number);
for i=1:2
    q_bar_i_fission=q_bar_i_fission+chi(:,:,g).*nu_sigma_f(:,:,i).*flux_old(:,:,i)/k_old;
end
q_bar_i_scattering_ind=zeros(mesh_center_abscissa_number,mesh_center_ordinate_number);
for i=1:2
    if(i~=g)
        q_bar_i_scattering_ind=q_bar_i_scattering_ind+sigma_s(:,:,i,g).*flux_old(:,:,i);
    end
end
q_bar_i_scattering_self=sigma_s(:,:,g,g).*flux_old(:,:,g);

q_bar_i=q_bar_i_fission+q_bar_i_scattering_ind+q_bar_i_scattering_self;

%%
Q_i_x_fission=zeros(mesh_center_abscissa_number,mesh_center_ordinate_number);
for i=1:2
    Q_i_x_fission=Q_i_x_fission+chi(:,:,g).*nu_sigma_f(:,:,i).*x_moment_flux(:,:,i)/k_old;
end

Q_i_x_scattering_ind=zeros(mesh_center_abscissa_number,mesh_center_ordinate_number);
for i=1:2
    if(i~=g)
        Q_i_x_scattering_ind=Q_i_x_scattering_ind+sigma_s(:,:,i,g).*x_moment_flux(:,:,i);
    end
end
Q_i_x_scattering_self=sigma_s(:,:,g,g).*x_moment_flux(:,:,g);

Q_i_x=Q_i_x_fission+Q_i_x_scattering_ind+Q_i_x_scattering_self;

%%
Q_i_y_fission=zeros(mesh_center_abscissa_number,mesh_center_ordinate_number);
for i=1:2
    Q_i_y_fission=Q_i_y_fission+chi(:,:,g).*nu_sigma_f(:,:,i).*y_moment_flux(:,:,i)/k_old;
end

Q_i_y_scattering_ind=zeros(mesh_center_abscissa_number,mesh_center_ordinate_number);
for i=1:2
    if(i~=g)
        Q_i_y_scattering_ind=Q_i_y_scattering_ind+sigma_s(:,:,i,g).*y_moment_flux(:,:,i);
    end
end
Q_i_y_scattering_self=sigma_s(:,:,g,g).*y_moment_flux(:,:,g);

Q_i_y=Q_i_y_fission+Q_i_y_scattering_ind+Q_i_y_scattering_self;



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
flux_new=zeros(mesh_center_abscissa_number,mesh_center_ordinate_number,2);

[flux_new(:,:,g),x_moment_flux(:,:,g),y_moment_flux(:,:,g),psi_bound]=transport_sweep(q_i_x,q_i_y,q_bar_i,c_i_xx,c_i_yy,c_i_xy,F_1,F_2,H,adj_len,s_len,ksi,alt_azim_theta,fin_d,x_c_t,y_c_t,X_i_c,Y_i_c,X,Y,dx,dy,N_a,sigma_t,psi_bound);

itrn=1;


while max(max(max(abs(flux_new-flux_old))))>tol
    flux_old=flux_new;

    q_bar_i_scattering_self=sigma_s(:,:,g,g).*flux_old(:,:,g);
    Q_i_x_scattering_self=sigma_s(:,:,g,g).*x_moment_flux(:,:,g);
    Q_i_y_scattering_self=sigma_s(:,:,g,g).*y_moment_flux(:,:,g);


    q_bar_i=q_bar_i_fission+q_bar_i_scattering_ind+q_bar_i_scattering_self;
    Q_i_x=Q_i_x_fission+Q_i_x_scattering_ind+Q_i_x_scattering_self;
    Q_i_y=Q_i_y_fission+Q_i_y_scattering_ind+Q_i_y_scattering_self;

    for j=1:mesh_center_ordinate_number
        for i=1:mesh_center_abscissa_number
            temp_M=M(:,:,i,j);
            solved=temp_M\[q_bar_i(i,j); Q_i_x(i,j); Q_i_y(i,j)];
            
            q_bar_i(i,j)=solved(1,1);
            q_i_x(i,j)=solved(2,1);
            q_i_y(i,j)=solved(3,1);
        end
    end

    [flux_new(:,:,g),x_moment_flux(:,:,g),y_moment_flux(:,:,g),psi_bound]=transport_sweep(q_i_x,q_i_y,q_bar_i,c_i_xx,c_i_yy,c_i_xy,F_1,F_2,H,adj_len,s_len,ksi,alt_azim_theta,fin_d,x_c_t,y_c_t,X_i_c,Y_i_c,X,Y,dx,dy,N_a,sigma_t,psi_bound);



    itrn=1+itrn
end

gth_flux=flux_new(:,:,g);
gth_x_moment_flux=x_moment_flux(:,:,g);
gth_y_moment_flux=y_moment_flux(:,:,g);