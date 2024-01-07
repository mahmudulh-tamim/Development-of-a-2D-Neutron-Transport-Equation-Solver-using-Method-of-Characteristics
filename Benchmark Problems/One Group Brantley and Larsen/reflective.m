function scaler_flux=reflective(S,exponential_portion,s_len,sum_s_len,alt_azim_theta,fin_d,X,Y,dx,dy,N_a,sigma_t,total_rays)
tol=10^(-7);

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


%%


weight_azimuthal=zeros(N_a,1);
weight_azimuthal(2:end-1,1)=0.5*(alt_azim_theta(3:end,1)-alt_azim_theta(1:end-2,1));
weight_azimuthal(1,1)=0.5*(alt_azim_theta(1,1)+alt_azim_theta(2,1));
weight_azimuthal(end,1)=0.5*(2*pi-alt_azim_theta(end-1,1)+alt_azim_theta(1,1));

%% end data


size_s_len=size(s_len);
psi_avg_old=zeros(size_s_len);

uuf=max(max(total_rays));
psi_bound=zeros(azimuthal_discretization_number,polar_discretization_number,uuf);

[psi_avg, psi_bound]=transport_sweep(S,exponential_portion,s_len,sum_s_len,alt_azim_theta,fin_d,X,Y,dx,dy,N_a,sigma_t,psi_bound);
it=1;
while max(max(max(max(abs(psi_avg-psi_avg_old)))))>tol
    psi_avg_old=psi_avg;
    [psi_avg, psi_bound]=transport_sweep(S,exponential_portion,s_len,sum_s_len,alt_azim_theta,fin_d,X,Y,dx,dy,N_a,sigma_t,psi_bound);
    it=it+1;
end
scaler_flux=zeros(mesh_center_ordinate_number,mesh_center_abscissa_number);

for j=1:mesh_center_ordinate_number
    for i=1:mesh_center_abscissa_number
        for az=1:azimuthal_discretization_number
            for p=1:polar_discretization_number
                scaler_flux(j,i)=psi_avg(j,i,az,p)*w(p,1)*weight_azimuthal(az,1)+scaler_flux(j,i);
            end
        end
    end
end
