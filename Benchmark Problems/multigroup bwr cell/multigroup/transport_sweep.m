function [avg_psi,psi_bound]=transport_sweep(S,exponential_portion,s_len,sum_s_len,alt_azim_theta,fin_d,X,Y,dx,dy,N_a,sigma_t,psi_bound)





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

ray_index_count=zeros(mesh_center_ordinate_number,mesh_center_abscissa_number,azimuthal_discretization_number,polar_discretization_number);
total_rays=zeros(azimuthal_discretization_number,polar_discretization_number);

psi_in=zeros(mesh_center_ordinate_number,mesh_center_abscissa_number,azimuthal_discretization_number,polar_discretization_number,100);
psi_out=zeros(mesh_center_ordinate_number,mesh_center_abscissa_number,azimuthal_discretization_number,polar_discretization_number,100);
del_psi=zeros(mesh_center_ordinate_number,mesh_center_abscissa_number,azimuthal_discretization_number,polar_discretization_number,100);
avg_psi=zeros(mesh_center_ordinate_number,mesh_center_abscissa_number,azimuthal_discretization_number,polar_discretization_number);

scaler_flux=zeros(mesh_center_ordinate_number,mesh_center_abscissa_number);



%% bottom to top rays


%left to right angle less than pi/2
for az_count=1:N_a/4
    
    ray_spacing_x=fin_d(az_count,1)*abs(csc(alt_azim_theta(az_count,1)));
    ray_spacing_y=fin_d(az_count,1)*abs(sec(alt_azim_theta(az_count,1)));
    ray_pos_x_bound=(ray_spacing_x/2:ray_spacing_x:X)';
    ray_pos_y_bound=(ray_spacing_y/2:ray_spacing_y:Y)';


    num_x_rays=length(ray_pos_x_bound);

    num_y_rays=length(ray_pos_y_bound);

    
    for pol_count=1:polar_discretization_number
        total_rays(az_count,pol_count)=num_y_rays+num_x_rays;
        ray_iden=total_rays(az_count,pol_count);
            i_x=1;
        
            for p_y=num_y_rays:-1:1
            
                if abs(ceil(ray_pos_y_bound(p_y,1)/dy)-ray_pos_y_bound(p_y,1)/dy)>10^(-14)
                    i_y=ceil(ray_pos_y_bound(p_y,1)/dy);
                else
                    i_y=ceil(ray_pos_y_bound(p_y,1)/dy)+1;
                end
                
                in_dx=i_x;
                in_dy=i_y;
        
                x_old=0;
                y_old=ray_pos_y_bound(p_y,1);
        
                while in_dx<=mesh_center_abscissa_number && in_dy<=mesh_center_ordinate_number
                    
                    ray_index_count(in_dy,in_dx,az_count,pol_count)=ray_index_count(in_dy,in_dx,az_count,pol_count)+1;


                    if in_dx==1
                          psi_in(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count))=psi_bound(az_count, pol_count, ray_iden);
                    else
                          psi_in(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count))=psi_out(in_dy_old,in_dx_old,az_count,pol_count,ray_index_count(in_dy_old,in_dx_old,az_count,pol_count));
    
                    end
                    del_psi(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count))=(psi_in(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count))-S(in_dy,in_dx)/sigma_t(in_dy, in_dx))*exponential_portion(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count));
                    psi_out(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count))=psi_in(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count))-del_psi(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count));
                    avg_psi(in_dy,in_dx,az_count,pol_count)=avg_psi(in_dy,in_dx,az_count,pol_count)+del_psi(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count))/(sigma_t(in_dy, in_dx)*sum_s_len(in_dy,in_dx,az_count,pol_count))+S(in_dy,in_dx)*s_len(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count))/(sigma_t(in_dy, in_dx)*sum_s_len(in_dy,in_dx,az_count,pol_count));
                    scaler_flux(in_dy,in_dx)=scaler_flux(in_dy,in_dx)+(del_psi(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count))/(sigma_t(in_dy, in_dx)*sum_s_len(in_dy,in_dx,az_count,pol_count))+S(in_dy,in_dx)*s_len(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count))/(sigma_t(in_dy, in_dx)*sum_s_len(in_dy,in_dx,az_count,pol_count)))*weight_azimuthal(az_count,1)*w(pol_count,1);

                    
                    x_new=dx*in_dx;
                    y_new=tan(alt_azim_theta(az_count,1))*(x_new-x_old)+y_old;
        
                    if y_new>dy*in_dy && abs(y_new-dy*in_dy)>10^(-14)
                        y_new=dy*in_dy;
                        x_new=x_old+(y_new-y_old)/tan(alt_azim_theta(az_count,1));
    
                        
                        
                        in_dx_old=in_dx;
                        in_dy_old=in_dy;

                        in_dy=in_dy+1;
                        x_old=x_new;
                        y_old=y_new;
                    elseif abs(y_new-dy*in_dy)<=10^(-14)
    
                        
                        
                        in_dx_old=in_dx;
                        in_dy_old=in_dy;

                        in_dx=in_dx+1;
                        in_dy=in_dy+1;
                        x_old=x_new;
                        y_old=y_new; 
                    else
                        
                        in_dx_old=in_dx;
                        in_dy_old=in_dy;

                        in_dx=in_dx+1;
                        x_old=x_new;
                        y_old=y_new;
                    end

                    if(in_dy==mesh_center_ordinate_number+1)
                       psi_bound(azimuthal_discretization_number-az_count+1,pol_count,ray_iden-num_y_rays)=psi_out(in_dy_old,in_dx_old,az_count,pol_count,ray_index_count(in_dy_old,in_dx_old,az_count,pol_count));
                    elseif in_dx==mesh_center_abscissa_number+1
                       psi_bound(azimuthal_discretization_number/2-az_count+1,pol_count,ray_iden+num_x_rays)=psi_out(in_dy_old,in_dx_old,az_count,pol_count,ray_index_count(in_dy_old,in_dx_old,az_count,pol_count));

                    end
                end
                ray_iden=ray_iden-1;
            end
        
            i_y=1;
        
            for p_x=1:num_x_rays
        
                if abs(ceil(ray_pos_x_bound(p_x,1)/dx)-ray_pos_x_bound(p_x,1)/dx)>10^(-14)
                    i_x=ceil(ray_pos_x_bound(p_x,1)/dx);
                else
                    i_x=ceil(ray_pos_x_bound(p_x,1)/dx)+1;
                end
                
                in_dx=i_x;
                in_dy=i_y;
        
                x_old=ray_pos_x_bound(p_x,1);
                y_old=0;
        
                while in_dx<=mesh_center_abscissa_number && in_dy<=mesh_center_ordinate_number
                    ray_index_count(in_dy,in_dx,az_count,pol_count)=ray_index_count(in_dy,in_dx,az_count,pol_count)+1;
                    if in_dy==1
                            psi_in(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count))=psi_bound(az_count, pol_count, ray_iden);
                    else
                            psi_in(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count))=psi_out(in_dy_old,in_dx_old,az_count,pol_count,ray_index_count(in_dy_old,in_dx_old,az_count,pol_count));
    
                    end
    
                    del_psi(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count))=(psi_in(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count))-S(in_dy,in_dx)/sigma_t(in_dy, in_dx))*exponential_portion(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count));
                    psi_out(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count))=psi_in(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count))-del_psi(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count));
    
                   avg_psi(in_dy,in_dx,az_count,pol_count)=avg_psi(in_dy,in_dx,az_count,pol_count)+del_psi(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count))/(sigma_t(in_dy, in_dx)*sum_s_len(in_dy,in_dx,az_count,pol_count))+S(in_dy,in_dx)*s_len(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count))/(sigma_t(in_dy, in_dx)*sum_s_len(in_dy,in_dx,az_count,pol_count));

    
                   scaler_flux(in_dy,in_dx)=scaler_flux(in_dy,in_dx)+(del_psi(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count))/(sigma_t(in_dy, in_dx)*sum_s_len(in_dy,in_dx,az_count,pol_count))+S(in_dy,in_dx)*s_len(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count))/(sigma_t(in_dy, in_dx)*sum_s_len(in_dy,in_dx,az_count,pol_count)))*weight_azimuthal(az_count,1)*w(pol_count,1);
                        
                    x_new=dx*in_dx;
                    y_new=tan(alt_azim_theta(az_count,1))*(x_new-x_old)+y_old;
        
                    if y_new>dy*in_dy && abs(y_new-dy*in_dy)>10^(-14)
                        y_new=dy*in_dy;
                        x_new=x_old+(y_new-y_old)/tan(alt_azim_theta(az_count,1));
                        
                        
                        in_dx_old=in_dx;
                        in_dy_old=in_dy;

                        in_dy=in_dy+1;
                        x_old=x_new;
                        y_old=y_new;
                    elseif abs(y_new-dy*in_dy)<=10^(-14)
                        
                        
                        in_dx_old=in_dx;
                        in_dy_old=in_dy;

                        in_dx=in_dx+1;
                        in_dy=in_dy+1;
                        x_old=x_new;
                        y_old=y_new; 
                     
                    else
                        
                        in_dx_old=in_dx;
                        in_dy_old=in_dy;

                        in_dx=in_dx+1;
                        x_old=x_new;
                        y_old=y_new;
                    end
                    if(in_dy==mesh_center_ordinate_number+1)
                       psi_bound(azimuthal_discretization_number-az_count+1,pol_count,ray_iden-num_y_rays)=psi_out(in_dy_old,in_dx_old,az_count,pol_count,ray_index_count(in_dy_old,in_dx_old,az_count,pol_count));
                    elseif in_dx==mesh_center_abscissa_number+1
                       psi_bound(azimuthal_discretization_number/2-az_count+1,pol_count,ray_iden+num_x_rays)=psi_out(in_dy_old,in_dx_old,az_count,pol_count,ray_index_count(in_dy_old,in_dx_old,az_count,pol_count));

                    end
                end
                ray_iden=ray_iden-1;
            end
    end
end


     
%right to left angle greater than pi/2 less than pi
for az_count=N_a/4+1:N_a/2
   
    ray_spacing_x=fin_d(az_count,1)*abs(csc(alt_azim_theta(az_count,1)));
    ray_spacing_y=fin_d(az_count,1)*abs(sec(alt_azim_theta(az_count,1)));
    ray_pos_x_bound=(ray_spacing_x/2:ray_spacing_x:X)';
    ray_pos_y_bound=(ray_spacing_y/2:ray_spacing_y:Y)';

    num_x_rays=length(ray_pos_x_bound);
        
    num_y_rays=length(ray_pos_y_bound);
    

    for pol_count=1:polar_discretization_number
            
            total_rays(az_count,pol_count)=num_y_rays+num_x_rays;
            ray_iden=total_rays(az_count,pol_count);
        
            i_x=mesh_center_abscissa_number;
        
            for p_y=num_y_rays:-1:1
                
                if abs(ceil(ray_pos_y_bound(p_y,1)/dy)-ray_pos_y_bound(p_y,1)/dy)>10^(-14)
                    i_y=ceil(ray_pos_y_bound(p_y,1)/dy);
                else
                    i_y=ceil(ray_pos_y_bound(p_y,1)/dy)+1;
                end
        
                in_dx=i_x;
                in_dy=i_y;
        
                x_old=X;
                y_old=ray_pos_y_bound(p_y,1);
        
                while in_dx>=1 && in_dy<=mesh_center_ordinate_number
                     ray_index_count(in_dy,in_dx,az_count,pol_count)=ray_index_count(in_dy,in_dx,az_count,pol_count)+1;
                        if in_dx==mesh_center_abscissa_number
                            psi_in(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count))=psi_bound(az_count, pol_count, ray_iden);
                        else
                            psi_in(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count))=psi_out(in_dy_old,in_dx_old,az_count,pol_count,ray_index_count(in_dy_old,in_dx_old,az_count,pol_count));
    
                        end
    
                        del_psi(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count))=(psi_in(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count))-S(in_dy,in_dx)/sigma_t(in_dy, in_dx))*exponential_portion(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count));
                        psi_out(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count))=psi_in(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count))-del_psi(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count));
    
                        avg_psi(in_dy,in_dx,az_count,pol_count)=avg_psi(in_dy,in_dx,az_count,pol_count)+del_psi(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count))/(sigma_t(in_dy, in_dx)*sum_s_len(in_dy,in_dx,az_count,pol_count))+S(in_dy,in_dx)*s_len(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count))/(sigma_t(in_dy, in_dx)*sum_s_len(in_dy,in_dx,az_count,pol_count));

    
                        scaler_flux(in_dy,in_dx)=scaler_flux(in_dy,in_dx)+(del_psi(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count))/(sigma_t(in_dy, in_dx)*sum_s_len(in_dy,in_dx,az_count,pol_count))+S(in_dy,in_dx)*s_len(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count))/(sigma_t(in_dy, in_dx)*sum_s_len(in_dy,in_dx,az_count,pol_count)))*weight_azimuthal(az_count,1)*w(pol_count,1);
                        
                    x_new=dx*(in_dx-1);
                    y_new=tan(alt_azim_theta(az_count,1))*(x_new-x_old)+y_old;
        
                   if y_new>dy*in_dy && abs(y_new-dy*in_dy)>10^(-14)
                        y_new=dy*in_dy;
                        x_new=x_old+(y_new-y_old)/tan(alt_azim_theta(az_count,1));
                        
                        
                        in_dx_old=in_dx;
                        in_dy_old=in_dy;

                        in_dy=in_dy+1;
                        x_old=x_new;
                        y_old=y_new;
                    elseif abs(y_new-dy*in_dy)<=10^(-14)
                       
                        
                        in_dx_old=in_dx;
                        in_dy_old=in_dy;

                        in_dx=in_dx-1;
                        in_dy=in_dy+1;
                        x_old=x_new;
                        y_old=y_new; 

                    else
                        
                        
                        in_dx_old=in_dx;
                        in_dy_old=in_dy;

                        in_dx=in_dx-1;
                        x_old=x_new;
                        y_old=y_new;
                   end
                   if(in_dy==mesh_center_ordinate_number+1)
                       psi_bound(azimuthal_discretization_number-az_count+1,pol_count,ray_iden-num_y_rays)=psi_out(in_dy_old,in_dx_old,az_count,pol_count,ray_index_count(in_dy_old,in_dx_old,az_count,pol_count));
                    elseif in_dx==0
                       psi_bound(azimuthal_discretization_number/2-az_count+1,pol_count,ray_iden+num_x_rays)=psi_out(in_dy_old,in_dx_old,az_count,pol_count,ray_index_count(in_dy_old,in_dx_old,az_count,pol_count));

                    end
                end
                ray_iden=ray_iden-1;
            end
        
        
            i_y=1;
        
            for p_x=num_x_rays:-1:1
                if abs(floor(ray_pos_x_bound(p_x,1)/dx)-ray_pos_x_bound(p_x,1)/dx)>10^(-14)
                    i_x=ceil(ray_pos_x_bound(p_x,1)/dx);
                else 
                    i_x=floor(ray_pos_x_bound(p_x,1)/dx);
                end
        
                in_dx=i_x;
                in_dy=i_y;
        
                x_old=ray_pos_x_bound(p_x,1);
                y_old=0;
        
                while in_dx>=1 && in_dy<=mesh_center_ordinate_number
                    ray_index_count(in_dy,in_dx,az_count,pol_count)=ray_index_count(in_dy,in_dx,az_count,pol_count)+1;
                        if in_dy==1
                            psi_in(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count))=psi_bound(az_count, pol_count, ray_iden);
                        else
                            psi_in(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count))=psi_out(in_dy_old,in_dx_old,az_count,pol_count,ray_index_count(in_dy_old,in_dx_old,az_count,pol_count));
    
                        end
    
                        del_psi(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count))=(psi_in(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count))-S(in_dy,in_dx)/sigma_t(in_dy, in_dx))*exponential_portion(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count));
                        psi_out(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count))=psi_in(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count))-del_psi(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count));
    
                        avg_psi(in_dy,in_dx,az_count,pol_count)=avg_psi(in_dy,in_dx,az_count,pol_count)+del_psi(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count))/(sigma_t(in_dy, in_dx)*sum_s_len(in_dy,in_dx,az_count,pol_count))+S(in_dy,in_dx)*s_len(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count))/(sigma_t(in_dy, in_dx)*sum_s_len(in_dy,in_dx,az_count,pol_count));

    
                        scaler_flux(in_dy,in_dx)=scaler_flux(in_dy,in_dx)+(del_psi(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count))/(sigma_t(in_dy, in_dx)*sum_s_len(in_dy,in_dx,az_count,pol_count))+S(in_dy,in_dx)*s_len(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count))/(sigma_t(in_dy, in_dx)*sum_s_len(in_dy,in_dx,az_count,pol_count)))*weight_azimuthal(az_count,1)*w(pol_count,1);
                        

                    x_new=dx*(in_dx-1);
                  
                    y_new=tan(alt_azim_theta(az_count,1))*(x_new-x_old)+y_old;
        
                    if y_new>dy*in_dy && abs(y_new-dy*in_dy)>10^(-14)
                        y_new=dy*in_dy;
                        x_new=x_old+(y_new-y_old)/tan(alt_azim_theta(az_count,1));
                        
                        
                        in_dx_old=in_dx;
                        in_dy_old=in_dy;

                        in_dy=in_dy+1;
                        x_old=x_new;
                        y_old=y_new;
                    elseif abs(y_new-dy*in_dy)<=10^(-14)
                        
                        
                        in_dx_old=in_dx;
                        in_dy_old=in_dy;

                        in_dx=in_dx-1;
                        in_dy=in_dy+1;
                        x_old=x_new;
                        y_old=y_new; 

                    else
                        
                        
                        in_dx_old=in_dx;
                        in_dy_old=in_dy;

                        in_dx=in_dx-1;
                        x_old=x_new;
                        y_old=y_new;
                    end
                    if(in_dy==mesh_center_ordinate_number+1)
                       psi_bound(azimuthal_discretization_number-az_count+1,pol_count,ray_iden-num_y_rays)=psi_out(in_dy_old,in_dx_old,az_count,pol_count,ray_index_count(in_dy_old,in_dx_old,az_count,pol_count));
                    elseif in_dx==0
                       psi_bound(azimuthal_discretization_number/2-az_count+1,pol_count,ray_iden+num_x_rays)=psi_out(in_dy_old,in_dx_old,az_count,pol_count,ray_index_count(in_dy_old,in_dx_old,az_count,pol_count));

                    end
                end
                ray_iden=ray_iden-1;
            end
    end

end


%% top to bottom rays


%left to right angle less than 2*pi greater than 3*pi/2
for az_count=3*N_a/4+1:N_a
    
    ray_spacing_x=fin_d(az_count,1)*abs(csc(alt_azim_theta(az_count,1)));
    ray_spacing_y=fin_d(az_count,1)*abs(sec(alt_azim_theta(az_count,1)));
    ray_pos_x_bound=(ray_spacing_x/2:ray_spacing_x:X)';
    ray_pos_y_bound=(ray_spacing_y/2:ray_spacing_y:Y)';


    num_x_rays=length(ray_pos_x_bound);

    num_y_rays=length(ray_pos_y_bound);

    for pol_count=1:polar_discretization_number
            
            total_rays(az_count,pol_count)=num_y_rays+num_x_rays;
            ray_iden=total_rays(az_count,pol_count);

            i_x=1;
        
            for p_y=1:num_y_rays
        
                if abs(floor(ray_pos_y_bound(p_y,1)/dy)-ray_pos_y_bound(p_y,1)/dy)>10^(-14)
                    i_y=ceil(ray_pos_y_bound(p_y,1)/dy);
                else 
                    i_y=floor(ray_pos_y_bound(p_y,1)/dy);
                end
                
                in_dx=i_x;
                in_dy=i_y;
        
                x_old=0;
                y_old=ray_pos_y_bound(p_y,1);
        
                while in_dx<=mesh_center_abscissa_number && in_dy>=1
                    ray_index_count(in_dy,in_dx,az_count,pol_count)=ray_index_count(in_dy,in_dx,az_count,pol_count)+1;
                        if in_dx==1
                            psi_in(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count))=psi_bound(az_count, pol_count, ray_iden);
                        else
                            psi_in(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count))=psi_out(in_dy_old,in_dx_old,az_count,pol_count,ray_index_count(in_dy_old,in_dx_old,az_count,pol_count));
    
                        end
    
                        del_psi(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count))=(psi_in(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count))-S(in_dy,in_dx)/sigma_t(in_dy, in_dx))*exponential_portion(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count));
                        psi_out(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count))=psi_in(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count))-del_psi(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count));
    
                        avg_psi(in_dy,in_dx,az_count,pol_count)=avg_psi(in_dy,in_dx,az_count,pol_count)+del_psi(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count))/(sigma_t(in_dy, in_dx)*sum_s_len(in_dy,in_dx,az_count,pol_count))+S(in_dy,in_dx)*s_len(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count))/(sigma_t(in_dy, in_dx)*sum_s_len(in_dy,in_dx,az_count,pol_count));

    
                        scaler_flux(in_dy,in_dx)=scaler_flux(in_dy,in_dx)+(del_psi(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count))/(sigma_t(in_dy, in_dx)*sum_s_len(in_dy,in_dx,az_count,pol_count))+S(in_dy,in_dx)*s_len(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count))/(sigma_t(in_dy, in_dx)*sum_s_len(in_dy,in_dx,az_count,pol_count)))*weight_azimuthal(az_count,1)*w(pol_count,1);
                        
                    x_new=dx*in_dx;
                    y_new=tan(alt_azim_theta(az_count,1))*(x_new-x_old)+y_old;
        
                    if y_new<dy*(in_dy-1) && abs(y_new-dy*(in_dy-1))>10^(-14)
                        y_new=dy*(in_dy-1);
                        x_new=x_old+(y_new-y_old)/tan(alt_azim_theta(az_count,1));
                        
                        
                        in_dx_old=in_dx;
                        in_dy_old=in_dy;

                        in_dy=in_dy-1;
                        x_old=x_new;
                        y_old=y_new;
                    elseif abs(y_new-dy*(in_dy-1))<=10^(-14)
                        
                        
                        in_dx_old=in_dx;
                        in_dy_old=in_dy;

                        in_dx=in_dx+1;
                        in_dy=in_dy-1;
                        x_old=x_new;
                        y_old=y_new; 

                    else
                        
                        
                        in_dx_old=in_dx;
                        in_dy_old=in_dy;

                        in_dx=in_dx+1;
                        x_old=x_new;
                        y_old=y_new;
                    end
                    if(in_dy==0)
                       psi_bound(azimuthal_discretization_number-az_count+1,pol_count,ray_iden-num_y_rays)=psi_out(in_dy_old,in_dx_old,az_count,pol_count,ray_index_count(in_dy_old,in_dx_old,az_count,pol_count));
                    elseif in_dx==mesh_center_abscissa_number+1
                       psi_bound(3*azimuthal_discretization_number/2-az_count+1,pol_count,ray_iden+num_x_rays)=psi_out(in_dy_old,in_dx_old,az_count,pol_count,ray_index_count(in_dy_old,in_dx_old,az_count,pol_count));

                    end
                end
                ray_iden=ray_iden-1;
            end
        
            i_y=mesh_center_ordinate_number;
        
            for p_x=1:num_x_rays
        
                if abs(ceil(ray_pos_x_bound(p_x,1)/dx)-ray_pos_x_bound(p_x,1)/dx)>10^(-14)
                    i_x=ceil(ray_pos_x_bound(p_x,1)/dx);
                else
                    i_x=ceil(ray_pos_x_bound(p_x,1)/dx)+1;
                end
                
                in_dx=i_x;
                in_dy=i_y;
        
                x_old=ray_pos_x_bound(p_x,1);
                y_old=Y;
        
                while in_dx<=mesh_center_abscissa_number && in_dy>=1
                    ray_index_count(in_dy,in_dx,az_count,pol_count)=ray_index_count(in_dy,in_dx,az_count,pol_count)+1;
                        if in_dy==mesh_center_ordinate_number
                            psi_in(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count))=psi_bound(az_count, pol_count, ray_iden);
                        else
                            psi_in(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count))=psi_out(in_dy_old,in_dx_old,az_count,pol_count,ray_index_count(in_dy_old,in_dx_old,az_count,pol_count));
    
                        end
    
                        del_psi(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count))=(psi_in(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count))-S(in_dy,in_dx)/sigma_t(in_dy, in_dx))*exponential_portion(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count));
                        psi_out(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count))=psi_in(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count))-del_psi(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count));
    
                        avg_psi(in_dy,in_dx,az_count,pol_count)=avg_psi(in_dy,in_dx,az_count,pol_count)+del_psi(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count))/(sigma_t(in_dy, in_dx)*sum_s_len(in_dy,in_dx,az_count,pol_count))+S(in_dy,in_dx)*s_len(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count))/(sigma_t(in_dy, in_dx)*sum_s_len(in_dy,in_dx,az_count,pol_count));

    
                        scaler_flux(in_dy,in_dx)=scaler_flux(in_dy,in_dx)+(del_psi(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count))/(sigma_t(in_dy, in_dx)*sum_s_len(in_dy,in_dx,az_count,pol_count))+S(in_dy,in_dx)*s_len(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count))/(sigma_t(in_dy, in_dx)*sum_s_len(in_dy,in_dx,az_count,pol_count)))*weight_azimuthal(az_count,1)*w(pol_count,1);
                        

                    x_new=dx*in_dx;
                    y_new=tan(alt_azim_theta(az_count,1))*(x_new-x_old)+y_old;
        
                    if y_new<dy*(in_dy-1) && abs(y_new-dy*(in_dy-1))>10^(-14)
                        y_new=dy*(in_dy-1);
                        x_new=x_old+(y_new-y_old)/tan(alt_azim_theta(az_count,1));
                        
                        
                        in_dx_old=in_dx;
                        in_dy_old=in_dy;

                        in_dy=in_dy-1;
                        x_old=x_new;
                        y_old=y_new;
                    elseif abs(y_new-dy*(in_dy-1))<=10^(-14)
                        
                        in_dx_old=in_dx;
                        in_dy_old=in_dy;

                        in_dx=in_dx+1;
                        in_dy=in_dy-1;
                        x_old=x_new;
                        y_old=y_new; 
                    else
                        
                        in_dx_old=in_dx;
                        in_dy_old=in_dy;

                        in_dx=in_dx+1;
                        x_old=x_new;
                        y_old=y_new;
                    end
                    if(in_dy==0)
                       psi_bound(azimuthal_discretization_number-az_count+1,pol_count,ray_iden-num_y_rays)=psi_out(in_dy_old,in_dx_old,az_count,pol_count,ray_index_count(in_dy_old,in_dx_old,az_count,pol_count));
                    elseif in_dx==mesh_center_abscissa_number+1
                       psi_bound(3*azimuthal_discretization_number/2-az_count+1,pol_count,ray_iden+num_x_rays)=psi_out(in_dy_old,in_dx_old,az_count,pol_count,ray_index_count(in_dy_old,in_dx_old,az_count,pol_count));

                    end

                end
                ray_iden=ray_iden-1;
            end
    end
end

      
%right to left angle greater than pi less than 3pi/2
for az_count=N_a/2+1:3*N_a/4
    
    ray_spacing_x=fin_d(az_count,1)*abs(csc(alt_azim_theta(az_count,1)));
    ray_spacing_y=fin_d(az_count,1)*abs(sec(alt_azim_theta(az_count,1)));
    ray_pos_x_bound=(ray_spacing_x/2:ray_spacing_x:X)';
    ray_pos_y_bound=(ray_spacing_y/2:ray_spacing_y:Y)';


    num_x_rays=length(ray_pos_x_bound);

    num_y_rays=length(ray_pos_y_bound);
    for pol_count=1:polar_discretization_number
            
            total_rays(az_count,pol_count)=num_x_rays+num_y_rays;
            ray_iden=total_rays(az_count,pol_count);

            i_x=mesh_center_abscissa_number;
        
            for p_y=1:num_y_rays
                
                if abs(floor(ray_pos_y_bound(p_y,1)/dy)-ray_pos_y_bound(p_y,1)/dy)>10^(-14)
                    i_y=ceil(ray_pos_y_bound(p_y,1)/dy);
                else 
                    i_y=floor(ray_pos_y_bound(p_y,1)/dy);
                end
        
                in_dx=i_x;
                in_dy=i_y;
        
                x_old=X;
                y_old=ray_pos_y_bound(p_y,1);
        
                while in_dx>=1 && in_dy>=1
                    ray_index_count(in_dy,in_dx,az_count,pol_count)=ray_index_count(in_dy,in_dx,az_count,pol_count)+1;
                        if in_dx==mesh_center_abscissa_number
                            psi_in(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count))=psi_bound(az_count, pol_count, ray_iden);
                        else
                            psi_in(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count))=psi_out(in_dy_old,in_dx_old,az_count,pol_count,ray_index_count(in_dy_old,in_dx_old,az_count,pol_count));
    
                        end
    
                        del_psi(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count))=(psi_in(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count))-S(in_dy,in_dx)/sigma_t(in_dy, in_dx))*exponential_portion(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count));
                        psi_out(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count))=psi_in(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count))-del_psi(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count));
    
                        avg_psi(in_dy,in_dx,az_count,pol_count)=avg_psi(in_dy,in_dx,az_count,pol_count)+del_psi(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count))/(sigma_t(in_dy, in_dx)*sum_s_len(in_dy,in_dx,az_count,pol_count))+S(in_dy,in_dx)*s_len(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count))/(sigma_t(in_dy, in_dx)*sum_s_len(in_dy,in_dx,az_count,pol_count));

                        scaler_flux(in_dy,in_dx)=scaler_flux(in_dy,in_dx)+(del_psi(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count))/(sigma_t(in_dy, in_dx)*sum_s_len(in_dy,in_dx,az_count,pol_count))+S(in_dy,in_dx)*s_len(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count))/(sigma_t(in_dy, in_dx)*sum_s_len(in_dy,in_dx,az_count,pol_count)))*weight_azimuthal(az_count,1)*w(pol_count,1);
                        
                    x_new=dx*(in_dx-1);
                    y_new=tan(alt_azim_theta(az_count,1))*(x_new-x_old)+y_old;
        
                    if y_new<dy*(in_dy-1)&& abs(y_new-dy*(in_dy-1))>10^(-14)
                        y_new=dy*(in_dy-1);
                        x_new=x_old+(y_new-y_old)/tan(alt_azim_theta(az_count,1));
                         
    
                        
                        
                        
                        in_dx_old=in_dx;
                        in_dy_old=in_dy;

                        in_dy=in_dy-1;
                        x_old=x_new;
                        y_old=y_new;
                    elseif abs(y_new-dy*(in_dy-1))<=10^(-14)
                        
                        
                        in_dx_old=in_dx;
                        in_dy_old=in_dy;

                        in_dx=in_dx-1;
                        in_dy=in_dy-1;
                        x_old=x_new;
                        y_old=y_new; 

                    else
                        
                        
                        in_dx_old=in_dx;
                        in_dy_old=in_dy;

                        in_dx=in_dx-1;
                        x_old=x_new;
                        y_old=y_new;
                    end
                    if(in_dy==0)
                       psi_bound(azimuthal_discretization_number-az_count+1,pol_count,ray_iden-num_y_rays)=psi_out(in_dy_old,in_dx_old,az_count,pol_count,ray_index_count(in_dy_old,in_dx_old,az_count,pol_count));
                    elseif in_dx==0
                        
                       psi_bound(3*azimuthal_discretization_number/2-az_count+1,pol_count,ray_iden+num_x_rays)=psi_out(in_dy_old,in_dx_old,az_count,pol_count,ray_index_count(in_dy_old,in_dx_old,az_count,pol_count));

                    end
                end
                ray_iden=ray_iden-1;
            end
        
        
            i_y=mesh_center_ordinate_number;
        
            for p_x=num_x_rays:-1:1
                if abs(floor(ray_pos_x_bound(p_x,1)/dx)-ray_pos_x_bound(p_x,1)/dx)>10^(-14)
                    i_x=ceil(ray_pos_x_bound(p_x,1)/dx);
                else 
                    i_x=floor(ray_pos_x_bound(p_x,1)/dx);
                end
        
                in_dx=i_x;
                in_dy=i_y;
        
                x_old=ray_pos_x_bound(p_x,1);
                y_old=Y;
        
                while in_dx>=1 && in_dy>=1
                    ray_index_count(in_dy,in_dx,az_count,pol_count)=ray_index_count(in_dy,in_dx,az_count,pol_count)+1;
                        if in_dy==mesh_center_ordinate_number
                            
                            psi_in(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count))=psi_bound(az_count, pol_count, ray_iden);
                        else
                            psi_in(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count))=psi_out(in_dy_old,in_dx_old,az_count,pol_count,ray_index_count(in_dy_old,in_dx_old,az_count,pol_count));
    
                        end
    
                        del_psi(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count))=(psi_in(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count))-S(in_dy,in_dx)/sigma_t(in_dy, in_dx))*exponential_portion(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count));
                        psi_out(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count))=psi_in(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count))-del_psi(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count));
    
                        avg_psi(in_dy,in_dx,az_count,pol_count)=avg_psi(in_dy,in_dx,az_count,pol_count)+del_psi(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count))/(sigma_t(in_dy, in_dx)*sum_s_len(in_dy,in_dx,az_count,pol_count))+S(in_dy,in_dx)*s_len(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count))/(sigma_t(in_dy, in_dx)*sum_s_len(in_dy,in_dx,az_count,pol_count));

    
                        scaler_flux(in_dy,in_dx)=scaler_flux(in_dy,in_dx)+(del_psi(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count))/(sigma_t(in_dy, in_dx)*sum_s_len(in_dy,in_dx,az_count,pol_count))+S(in_dy,in_dx)*s_len(in_dy,in_dx,az_count,pol_count,ray_index_count(in_dy,in_dx,az_count,pol_count))/(sigma_t(in_dy, in_dx)*sum_s_len(in_dy,in_dx,az_count,pol_count)))*weight_azimuthal(az_count,1)*w(pol_count,1);
                        

                    x_new=dx*(in_dx-1);
                  
                    y_new=tan(alt_azim_theta(az_count,1))*(x_new-x_old)+y_old;
        
                    if y_new<dy*(in_dy-1) && abs(y_new-dy*(in_dy-1))>10^(-14)
                        y_new=dy*(in_dy-1);
                        x_new=x_old+(y_new-y_old)/tan(alt_azim_theta(az_count,1));
                        
                        
                        in_dx_old=in_dx;
                        in_dy_old=in_dy;

                        in_dy=in_dy-1;
                        x_old=x_new;
                        y_old=y_new;
                      elseif abs(y_new-dy*(in_dy-1))<=10^(-14)
                        
                        
                        in_dx_old=in_dx;
                        in_dy_old=in_dy;

                        in_dx=in_dx-1;
                        in_dy=in_dy-1;
                        x_old=x_new;
                        y_old=y_new; 

                    else
                        
                        
                        in_dx_old=in_dx;
                        in_dy_old=in_dy;

                        in_dx=in_dx-1;
                        x_old=x_new;
                        y_old=y_new;
                    end
                    if(in_dy==0)
                       psi_bound(azimuthal_discretization_number-az_count+1,pol_count,ray_iden-num_y_rays)=psi_out(in_dy_old,in_dx_old,az_count,pol_count,ray_index_count(in_dy_old,in_dx_old,az_count,pol_count));
                    elseif in_dx==0
                       psi_bound(3*azimuthal_discretization_number/2-az_count+1,pol_count,ray_iden+num_x_rays)=psi_out(in_dy_old,in_dx_old,az_count,pol_count,ray_index_count(in_dy_old,in_dx_old,az_count,pol_count));

                    end
                end
                ray_iden=ray_iden-1;
            end
    end

end



