function [F_1,F_2,H,ksi,x_c_t,y_c_t,X_i_c,Y_i_c,s_len,adj_len,alt_azim_theta,fin_d, mesh_center_abscissa_number,mesh_center_ordinate_number, total_rays,M,c_i_xx, c_i_yy,c_i_xy]=ray_tracing_for_other_data(X,Y,dx,dy,N_a,sigma_t)
%M_xx_i and M_yy_i and M_xy_i don't incorporate division by ksi or do they?


[F_1,F_2,G_2,H,ksi,x_c_t,y_c_t,X_i_c,Y_i_c,s_len,adj_len,alt_azim_theta,fin_d]=ray_tracing_polar_included(X,Y,dx,dy,N_a,sigma_t);









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

weight_azimuthal=zeros(N_a,1);
weight_azimuthal(2:end-1,1)=0.5*(alt_azim_theta(3:end,1)-alt_azim_theta(1:end-2,1));
weight_azimuthal(1,1)=0.5*(alt_azim_theta(1,1)+alt_azim_theta(2,1));
weight_azimuthal(end,1)=0.5*(2*pi-alt_azim_theta(end-1,1)+alt_azim_theta(1,1));


%initialization

ray_index_count=zeros(mesh_center_abscissa_number,mesh_center_ordinate_number,azimuthal_discretization_number,polar_discretization_number);




M_xx_i=zeros(mesh_center_abscissa_number,mesh_center_ordinate_number);
M_yy_i=zeros(mesh_center_abscissa_number,mesh_center_ordinate_number);
M_xy_i=zeros(mesh_center_abscissa_number,mesh_center_ordinate_number);

M=zeros(3,3,mesh_center_abscissa_number,mesh_center_ordinate_number);

c_i_xx=zeros(mesh_center_abscissa_number,mesh_center_ordinate_number);
c_i_yy=zeros(mesh_center_abscissa_number,mesh_center_ordinate_number);
c_i_xy=zeros(mesh_center_abscissa_number,mesh_center_ordinate_number);





total_rays=zeros(azimuthal_discretization_number,polar_discretization_number);


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
            i_x=1;

            total_rays(az_count,pol_count)=num_y_rays+num_x_rays;

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
                    ray_index_count(in_dx,in_dy,az_count,pol_count)=ray_index_count(in_dx,in_dy,az_count,pol_count)+1;
                   
                    omega_x=cos(alt_azim_theta(az_count,1))*mu(pol_count,1);
                    omega_y=sin(alt_azim_theta(az_count,1))*mu(pol_count,1);
                    x_c_track=x_c_t(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count,pol_count));
                    y_c_track=y_c_t(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count,pol_count));
                    
                    t_mki=s_len(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))/ksi(in_dx,in_dy,az_count,1);
                    M_xx_i(in_dx,in_dy)=M_xx_i(in_dx,in_dy)+1/(dx*dy)*mu(pol_count,1)*w(pol_count,1)*weight_azimuthal(az_count,1)*fin_d(az_count,1)*((omega_x)^2*t_mki^3/3+omega_x*(x_old-X_i_c(in_dx,in_dy))*t_mki^2+t_mki*(x_old-X_i_c(in_dx,in_dy))^2);
                    M_yy_i(in_dx,in_dy)=M_yy_i(in_dx,in_dy)+1/(dx*dy)*mu(pol_count,1)*w(pol_count,1)*weight_azimuthal(az_count,1)*fin_d(az_count,1)*((omega_y)^2*t_mki^3/3+omega_y*(y_old-Y_i_c(in_dx,in_dy))*t_mki^2+t_mki*(y_old-Y_i_c(in_dx,in_dy))^2);

                    M_xy_i(in_dx,in_dy)=M_xy_i(in_dx,in_dy)+1/(dx*dy)*mu(pol_count,1)*w(pol_count,1)*weight_azimuthal(az_count,1)*fin_d(az_count,1)*(omega_x*omega_y*t_mki^3/3+0.5*((x_old-X_i_c(in_dx,in_dy)*omega_y)+(y_old-Y_i_c(in_dx,in_dy))*omega_x)*t_mki^2+t_mki*(x_old-X_i_c(in_dx,in_dy))*(y_old-Y_i_c(in_dx,in_dy)));
                    
                    wp_bar=mu(pol_count,1)*w(pol_count,1);
                    wa_bar=weight_azimuthal(az_count,1)*fin_d(az_count,1);
                    cosphia=cos(alt_azim_theta(az_count,1));
                    sinphia=sin(alt_azim_theta(az_count,1));
                    temp_len=adj_len(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count,pol_count))/ksi(in_dx,in_dy,az_count);
                    temp_t_len=adj_len(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count,pol_count));
                    G_2_temp=G_2(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count));
                    c_i_xx(in_dx,in_dy)=c_i_xx(in_dx,in_dy)+1/(sigma_t(in_dx,in_dy)*dx*dy)*wp_bar*wa_bar*(cosphia*temp_len)^2*G_2_temp+1/(dx*dy)*w(pol_count,1)*wa_bar*(x_c_track-X_i_c(in_dx,in_dy))^2*temp_t_len;
                    c_i_yy(in_dx,in_dy)=c_i_yy(in_dx,in_dy)+1/(sigma_t(in_dx,in_dy)*dx*dy)*wp_bar*wa_bar*(sinphia*temp_len)^2*G_2_temp+1/(dx*dy)*w(pol_count,1)*wa_bar*(y_c_track-Y_i_c(in_dx,in_dy))^2*temp_t_len;

                    c_i_xy(in_dx,in_dy)=c_i_xy(in_dx,in_dy)+1/(sigma_t(in_dx,in_dy)*dx*dy)*wp_bar*wa_bar*sinphia*cosphia*(temp_len)^2*G_2_temp+1/(dx*dy)*wa_bar*w(pol_count,1)*(y_c_track-Y_i_c(in_dx,in_dy))*(x_c_track-X_i_c(in_dx,in_dy))*temp_t_len;

                    
                    x_new=dx*in_dx;
                    y_new=tan(alt_azim_theta(az_count,1))*(x_new-x_old)+y_old;

                    
        
                    if y_new>dy*in_dy && abs(y_new-dy*in_dy)>10^(-14)
                        y_new=dy*in_dy;
                        x_new=x_old+(y_new-y_old)/tan(alt_azim_theta(az_count,1));
                        
                        

                        in_dy=in_dy+1;
                        x_old=x_new;
                        y_old=y_new;
                    elseif abs(y_new-dy*in_dy)<=10^(-14)
                        
                        in_dx=in_dx+1;
                        in_dy=in_dy+1;
                        x_old=x_new;
                        y_old=y_new;
                    else
                        
                        in_dx=in_dx+1;
                        x_old=x_new;
                        y_old=y_new;
                    end
                end
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
                    ray_index_count(in_dx,in_dy,az_count,pol_count)=ray_index_count(in_dx,in_dy,az_count,pol_count)+1;
                    
                    omega_x=cos(alt_azim_theta(az_count,1))*mu(pol_count,1);
                    omega_y=sin(alt_azim_theta(az_count,1))*mu(pol_count,1);
                    x_c_track=x_c_t(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count,pol_count));
                    y_c_track=y_c_t(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count,pol_count));
                    t_mki=s_len(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))/ksi(in_dx,in_dy,az_count,1);
                    M_xx_i(in_dx,in_dy)=M_xx_i(in_dx,in_dy)+1/(dx*dy)*mu(pol_count,1)*w(pol_count,1)*weight_azimuthal(az_count,1)*fin_d(az_count,1)*((omega_x)^2*t_mki^3/3+omega_x*(x_old-X_i_c(in_dx,in_dy))*t_mki^2+t_mki*(x_old-X_i_c(in_dx,in_dy))^2);
                    M_yy_i(in_dx,in_dy)=M_yy_i(in_dx,in_dy)+1/(dx*dy)*mu(pol_count,1)*w(pol_count,1)*weight_azimuthal(az_count,1)*fin_d(az_count,1)*((omega_y)^2*t_mki^3/3+omega_y*(y_old-Y_i_c(in_dx,in_dy))*t_mki^2+t_mki*(y_old-Y_i_c(in_dx,in_dy))^2);

                    M_xy_i(in_dx,in_dy)=M_xy_i(in_dx,in_dy)+1/(dx*dy)*mu(pol_count,1)*w(pol_count,1)*weight_azimuthal(az_count,1)*fin_d(az_count,1)*(omega_x*omega_y*t_mki^3/3+0.5*((x_old-X_i_c(in_dx,in_dy)*omega_y)+(y_old-Y_i_c(in_dx,in_dy))*omega_x)*t_mki^2+t_mki*(x_old-X_i_c(in_dx,in_dy))*(y_old-Y_i_c(in_dx,in_dy)));
                    
                    wp_bar=mu(pol_count,1)*w(pol_count,1);
                    wa_bar=weight_azimuthal(az_count,1)*fin_d(az_count,1);
                    cosphia=cos(alt_azim_theta(az_count,1));
                    sinphia=sin(alt_azim_theta(az_count,1));
                    temp_len=adj_len(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count,pol_count))/ksi(in_dx,in_dy,az_count);
                    temp_t_len=adj_len(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count,pol_count));
                    G_2_temp=G_2(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count));
                    c_i_xx(in_dx,in_dy)=c_i_xx(in_dx,in_dy)+1/(sigma_t(in_dx,in_dy)*dx*dy)*wp_bar*wa_bar*(cosphia*temp_len)^2*G_2_temp+1/(dx*dy)*w(pol_count,1)*wa_bar*(x_c_track-X_i_c(in_dx,in_dy))^2*temp_t_len;
                    c_i_yy(in_dx,in_dy)=c_i_yy(in_dx,in_dy)+1/(sigma_t(in_dx,in_dy)*dx*dy)*wp_bar*wa_bar*(sinphia*temp_len)^2*G_2_temp+1/(dx*dy)*w(pol_count,1)*wa_bar*(y_c_track-Y_i_c(in_dx,in_dy))^2*temp_t_len;

                    c_i_xy(in_dx,in_dy)=c_i_xy(in_dx,in_dy)+1/(sigma_t(in_dx,in_dy)*dx*dy)*wp_bar*wa_bar*sinphia*cosphia*(temp_len)^2*G_2_temp+1/(dx*dy)*wa_bar*w(pol_count,1)*(y_c_track-Y_i_c(in_dx,in_dy))*(x_c_track-X_i_c(in_dx,in_dy))*temp_t_len;

                    
                    x_new=dx*in_dx;
                    y_new=tan(alt_azim_theta(az_count,1))*(x_new-x_old)+y_old;
        
                    if y_new>dy*in_dy && abs(y_new-dy*in_dy)>10^(-14)
                        y_new=dy*in_dy;
                        x_new=x_old+(y_new-y_old)/tan(alt_azim_theta(az_count,1));

                        
                        in_dy=in_dy+1;
                        x_old=x_new;
                        y_old=y_new;
                     elseif abs(y_new-dy*in_dy)<=10^(-14)
                        
                        in_dx=in_dx+1;
                        in_dy=in_dy+1;
                        x_old=x_new;
                        y_old=y_new;

                    else
                        
                        in_dx=in_dx+1;
                        x_old=x_new;
                        y_old=y_new;
                    end
                end
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
                    ray_index_count(in_dx,in_dy,az_count,pol_count)=ray_index_count(in_dx,in_dy,az_count,pol_count)+1;
                    omega_x=cos(alt_azim_theta(az_count,1))*mu(pol_count,1);
                    omega_y=sin(alt_azim_theta(az_count,1))*mu(pol_count,1);
                    x_c_track=x_c_t(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count,pol_count));
                    y_c_track=y_c_t(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count,pol_count));
                    t_mki=s_len(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))/ksi(in_dx,in_dy,az_count,1);
                    M_xx_i(in_dx,in_dy)=M_xx_i(in_dx,in_dy)+1/(dx*dy)*mu(pol_count,1)*w(pol_count,1)*weight_azimuthal(az_count,1)*fin_d(az_count,1)*((omega_x)^2*t_mki^3/3+omega_x*(x_old-X_i_c(in_dx,in_dy))*t_mki^2+t_mki*(x_old-X_i_c(in_dx,in_dy))^2);
                    M_yy_i(in_dx,in_dy)=M_yy_i(in_dx,in_dy)+1/(dx*dy)*mu(pol_count,1)*w(pol_count,1)*weight_azimuthal(az_count,1)*fin_d(az_count,1)*((omega_y)^2*t_mki^3/3+omega_y*(y_old-Y_i_c(in_dx,in_dy))*t_mki^2+t_mki*(y_old-Y_i_c(in_dx,in_dy))^2);

                    M_xy_i(in_dx,in_dy)=M_xy_i(in_dx,in_dy)+1/(dx*dy)*mu(pol_count,1)*w(pol_count,1)*weight_azimuthal(az_count,1)*fin_d(az_count,1)*(omega_x*omega_y*t_mki^3/3+0.5*((x_old-X_i_c(in_dx,in_dy)*omega_y)+(y_old-Y_i_c(in_dx,in_dy))*omega_x)*t_mki^2+t_mki*(x_old-X_i_c(in_dx,in_dy))*(y_old-Y_i_c(in_dx,in_dy)));

                   
                    x_new=dx*(in_dx-1);
                    y_new=tan(alt_azim_theta(az_count,1))*(x_new-x_old)+y_old;
        
                    if y_new>dy*in_dy && abs(y_new-dy*in_dy)>10^(-14)
                        y_new=dy*in_dy;
                        x_new=x_old+(y_new-y_old)/tan(alt_azim_theta(az_count,1));
                         
                        in_dy=in_dy+1;
                        x_old=x_new;
                        y_old=y_new;
                     elseif abs(y_new-dy*in_dy)<=10^(-14)
                       
                        in_dx=in_dx-1;
                        in_dy=in_dy+1;
                        x_old=x_new;
                        y_old=y_new;

                    else
                        
                        in_dx=in_dx-1;
                        x_old=x_new;
                        y_old=y_new;
                    end
                end
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
                    ray_index_count(in_dx,in_dy,az_count,pol_count)=ray_index_count(in_dx,in_dy,az_count,pol_count)+1;
                    omega_x=cos(alt_azim_theta(az_count,1))*mu(pol_count,1);
                    omega_y=sin(alt_azim_theta(az_count,1))*mu(pol_count,1);
                    x_c_track=x_c_t(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count,pol_count));
                    y_c_track=y_c_t(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count,pol_count));
                    t_mki=s_len(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))/ksi(in_dx,in_dy,az_count,1);
                    M_xx_i(in_dx,in_dy)=M_xx_i(in_dx,in_dy)+1/(dx*dy)*mu(pol_count,1)*w(pol_count,1)*weight_azimuthal(az_count,1)*fin_d(az_count,1)*((omega_x)^2*t_mki^3/3+omega_x*(x_old-X_i_c(in_dx,in_dy))*t_mki^2+t_mki*(x_old-X_i_c(in_dx,in_dy))^2);
                    M_yy_i(in_dx,in_dy)=M_yy_i(in_dx,in_dy)+1/(dx*dy)*mu(pol_count,1)*w(pol_count,1)*weight_azimuthal(az_count,1)*fin_d(az_count,1)*((omega_y)^2*t_mki^3/3+omega_y*(y_old-Y_i_c(in_dx,in_dy))*t_mki^2+t_mki*(y_old-Y_i_c(in_dx,in_dy))^2);

                    M_xy_i(in_dx,in_dy)=M_xy_i(in_dx,in_dy)+1/(dx*dy)*mu(pol_count,1)*w(pol_count,1)*weight_azimuthal(az_count,1)*fin_d(az_count,1)*(omega_x*omega_y*t_mki^3/3+0.5*((x_old-X_i_c(in_dx,in_dy)*omega_y)+(y_old-Y_i_c(in_dx,in_dy))*omega_x)*t_mki^2+t_mki*(x_old-X_i_c(in_dx,in_dy))*(y_old-Y_i_c(in_dx,in_dy)));
                    
                    

                    x_new=dx*(in_dx-1);
                  
                    y_new=tan(alt_azim_theta(az_count,1))*(x_new-x_old)+y_old;
        
                    if y_new>dy*in_dy && abs(y_new-dy*in_dy)>10^(-14)
                        y_new=dy*in_dy;
                        x_new=x_old+(y_new-y_old)/tan(alt_azim_theta(az_count,1));
                        
                        in_dy=in_dy+1;
                        x_old=x_new;
                        y_old=y_new;
                     elseif abs(y_new-dy*in_dy)<=10^(-14)
                        
                        in_dx=in_dx-1;
                        in_dy=in_dy+1;
                        x_old=x_new;
                        y_old=y_new;

                    else
                        
                        in_dx=in_dx-1;
                        x_old=x_new;
                        y_old=y_new;
                    end
                end
            end
    end

end


%% top to bottom rays


%left to right angle less than pi/2
for az_count=3*N_a/4+1:N_a
    
    ray_spacing_x=fin_d(az_count,1)*abs(csc(alt_azim_theta(az_count,1)));
    ray_spacing_y=fin_d(az_count,1)*abs(sec(alt_azim_theta(az_count,1)));
    ray_pos_x_bound=(ray_spacing_x/2:ray_spacing_x:X)';
    ray_pos_y_bound=(ray_spacing_y/2:ray_spacing_y:Y)';


    num_x_rays=length(ray_pos_x_bound);

    num_y_rays=length(ray_pos_y_bound);

    for pol_count=1:polar_discretization_number
            total_rays(az_count,pol_count)=num_y_rays+num_x_rays;
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
                    ray_index_count(in_dx,in_dy,az_count,pol_count)=ray_index_count(in_dx,in_dy,az_count,pol_count)+1;
                    omega_x=cos(alt_azim_theta(az_count,1))*mu(pol_count,1);
                    omega_y=sin(alt_azim_theta(az_count,1))*mu(pol_count,1);
                    x_c_track=x_c_t(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count,pol_count));
                    y_c_track=y_c_t(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count,pol_count));
                    t_mki=s_len(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))/ksi(in_dx,in_dy,az_count,1);
                    M_xx_i(in_dx,in_dy)=M_xx_i(in_dx,in_dy)+1/(dx*dy)*mu(pol_count,1)*w(pol_count,1)*weight_azimuthal(az_count,1)*fin_d(az_count,1)*((omega_x)^2*t_mki^3/3+omega_x*(x_old-X_i_c(in_dx,in_dy))*t_mki^2+t_mki*(x_old-X_i_c(in_dx,in_dy))^2);
                    M_yy_i(in_dx,in_dy)=M_yy_i(in_dx,in_dy)+1/(dx*dy)*mu(pol_count,1)*w(pol_count,1)*weight_azimuthal(az_count,1)*fin_d(az_count,1)*((omega_y)^2*t_mki^3/3+omega_y*(y_old-Y_i_c(in_dx,in_dy))*t_mki^2+t_mki*(y_old-Y_i_c(in_dx,in_dy))^2);

                    M_xy_i(in_dx,in_dy)=M_xy_i(in_dx,in_dy)+1/(dx*dy)*mu(pol_count,1)*w(pol_count,1)*weight_azimuthal(az_count,1)*fin_d(az_count,1)*(omega_x*omega_y*t_mki^3/3+0.5*((x_old-X_i_c(in_dx,in_dy)*omega_y)+(y_old-Y_i_c(in_dx,in_dy))*omega_x)*t_mki^2+t_mki*(x_old-X_i_c(in_dx,in_dy))*(y_old-Y_i_c(in_dx,in_dy)));
                   
                    wp_bar=mu(pol_count,1)*w(pol_count,1);
                    wa_bar=weight_azimuthal(az_count,1)*fin_d(az_count,1);
                    cosphia=cos(alt_azim_theta(az_count,1));
                    sinphia=sin(alt_azim_theta(az_count,1));
                    temp_len=adj_len(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count,pol_count))/ksi(in_dx,in_dy,az_count);
                    temp_t_len=adj_len(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count,pol_count));
                    G_2_temp=G_2(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count));
                    c_i_xx(in_dx,in_dy)=c_i_xx(in_dx,in_dy)+1/(sigma_t(in_dx,in_dy)*dx*dy)*wp_bar*wa_bar*(cosphia*temp_len)^2*G_2_temp+1/(dx*dy)*w(pol_count,1)*wa_bar*(x_c_track-X_i_c(in_dx,in_dy))^2*temp_t_len;
                    c_i_yy(in_dx,in_dy)=c_i_yy(in_dx,in_dy)+1/(sigma_t(in_dx,in_dy)*dx*dy)*wp_bar*wa_bar*(sinphia*temp_len)^2*G_2_temp+1/(dx*dy)*w(pol_count,1)*wa_bar*(y_c_track-Y_i_c(in_dx,in_dy))^2*temp_t_len;

                    c_i_xy(in_dx,in_dy)=c_i_xy(in_dx,in_dy)+1/(sigma_t(in_dx,in_dy)*dx*dy)*wp_bar*wa_bar*sinphia*cosphia*(temp_len)^2*G_2_temp+1/(dx*dy)*wa_bar*w(pol_count,1)*(y_c_track-Y_i_c(in_dx,in_dy))*(x_c_track-X_i_c(in_dx,in_dy))*temp_t_len;

                    

                    x_new=dx*in_dx;
                    y_new=tan(alt_azim_theta(az_count,1))*(x_new-x_old)+y_old;
        
                    if y_new<dy*(in_dy-1) && abs(y_new-dy*(in_dy-1))>10^(-14)
                        y_new=dy*(in_dy-1);
                        x_new=x_old+(y_new-y_old)/tan(alt_azim_theta(az_count,1));
                        
                        in_dy=in_dy-1;
                        x_old=x_new;
                        y_old=y_new;

                     elseif abs(y_new-dy*(in_dy-1))<=10^(-14)
                        
                        in_dx=in_dx+1;
                        in_dy=in_dy-1;
                        x_old=x_new;
                        y_old=y_new;
                    else
                        
                        in_dx=in_dx+1;
                        x_old=x_new;
                        y_old=y_new;
                    end
                end
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
                    ray_index_count(in_dx,in_dy,az_count,pol_count)=ray_index_count(in_dx,in_dy,az_count,pol_count)+1;
                    omega_x=cos(alt_azim_theta(az_count,1))*mu(pol_count,1);
                    omega_y=sin(alt_azim_theta(az_count,1))*mu(pol_count,1);
                    x_c_track=x_c_t(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count,pol_count));
                    y_c_track=y_c_t(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count,pol_count));
                    t_mki=s_len(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))/ksi(in_dx,in_dy,az_count,1);
                    M_xx_i(in_dx,in_dy)=M_xx_i(in_dx,in_dy)+1/(dx*dy)*mu(pol_count,1)*w(pol_count,1)*weight_azimuthal(az_count,1)*fin_d(az_count,1)*((omega_x)^2*t_mki^3/3+omega_x*(x_old-X_i_c(in_dx,in_dy))*t_mki^2+t_mki*(x_old-X_i_c(in_dx,in_dy))^2);
                    M_yy_i(in_dx,in_dy)=M_yy_i(in_dx,in_dy)+1/(dx*dy)*mu(pol_count,1)*w(pol_count,1)*weight_azimuthal(az_count,1)*fin_d(az_count,1)*((omega_y)^2*t_mki^3/3+omega_y*(y_old-Y_i_c(in_dx,in_dy))*t_mki^2+t_mki*(y_old-Y_i_c(in_dx,in_dy))^2);

                    M_xy_i(in_dx,in_dy)=M_xy_i(in_dx,in_dy)+1/(dx*dy)*mu(pol_count,1)*w(pol_count,1)*weight_azimuthal(az_count,1)*fin_d(az_count,1)*(omega_x*omega_y*t_mki^3/3+0.5*((x_old-X_i_c(in_dx,in_dy)*omega_y)+(y_old-Y_i_c(in_dx,in_dy))*omega_x)*t_mki^2+t_mki*(x_old-X_i_c(in_dx,in_dy))*(y_old-Y_i_c(in_dx,in_dy)));
                    
                    wp_bar=mu(pol_count,1)*w(pol_count,1);
                    wa_bar=weight_azimuthal(az_count,1)*fin_d(az_count,1);
                    cosphia=cos(alt_azim_theta(az_count,1));
                    sinphia=sin(alt_azim_theta(az_count,1));
                    temp_len=adj_len(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count,pol_count))/ksi(in_dx,in_dy,az_count);
                    temp_t_len=adj_len(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count,pol_count));
                    G_2_temp=G_2(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count));
                    c_i_xx(in_dx,in_dy)=c_i_xx(in_dx,in_dy)+1/(sigma_t(in_dx,in_dy)*dx*dy)*wp_bar*wa_bar*(cosphia*temp_len)^2*G_2_temp+1/(dx*dy)*w(pol_count,1)*wa_bar*(x_c_track-X_i_c(in_dx,in_dy))^2*temp_t_len;
                    c_i_yy(in_dx,in_dy)=c_i_yy(in_dx,in_dy)+1/(sigma_t(in_dx,in_dy)*dx*dy)*wp_bar*wa_bar*(sinphia*temp_len)^2*G_2_temp+1/(dx*dy)*w(pol_count,1)*wa_bar*(y_c_track-Y_i_c(in_dx,in_dy))^2*temp_t_len;

                    c_i_xy(in_dx,in_dy)=c_i_xy(in_dx,in_dy)+1/(sigma_t(in_dx,in_dy)*dx*dy)*wp_bar*wa_bar*sinphia*cosphia*(temp_len)^2*G_2_temp+1/(dx*dy)*wa_bar*w(pol_count,1)*(y_c_track-Y_i_c(in_dx,in_dy))*(x_c_track-X_i_c(in_dx,in_dy))*temp_t_len;

                    

                    x_new=dx*in_dx;
                    y_new=tan(alt_azim_theta(az_count,1))*(x_new-x_old)+y_old;
        
                    if y_new<dy*(in_dy-1) && abs(y_new-dy*(in_dy-1))>10^(-14)
                        y_new=dy*(in_dy-1);
                        x_new=x_old+(y_new-y_old)/tan(alt_azim_theta(az_count,1));
                        
                        in_dy=in_dy-1;
                        x_old=x_new;
                        y_old=y_new;
                     elseif abs(y_new-dy*(in_dy-1))<=10^(-14)
                        
                        in_dx=in_dx+1;
                        in_dy=in_dy-1;
                        x_old=x_new;
                        y_old=y_new;
                    else
                        
                        in_dx=in_dx+1;
                        x_old=x_new;
                        y_old=y_new;
                    end
                end
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
            total_rays(az_count,pol_count)=num_y_rays+num_x_rays;
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
                    ray_index_count(in_dx,in_dy,az_count,pol_count)=ray_index_count(in_dx,in_dy,az_count,pol_count)+1;
                    omega_x=cos(alt_azim_theta(az_count,1))*mu(pol_count,1);
                    omega_y=sin(alt_azim_theta(az_count,1))*mu(pol_count,1);
                    x_c_track=x_c_t(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count,pol_count));
                    y_c_track=y_c_t(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count,pol_count));
                    t_mki=s_len(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))/ksi(in_dx,in_dy,az_count,1);
                    M_xx_i(in_dx,in_dy)=M_xx_i(in_dx,in_dy)+1/(dx*dy)*mu(pol_count,1)*w(pol_count,1)*weight_azimuthal(az_count,1)*fin_d(az_count,1)*((omega_x)^2*t_mki^3/3+omega_x*(x_old-X_i_c(in_dx,in_dy))*t_mki^2+t_mki*(x_old-X_i_c(in_dx,in_dy))^2);
                    M_yy_i(in_dx,in_dy)=M_yy_i(in_dx,in_dy)+1/(dx*dy)*mu(pol_count,1)*w(pol_count,1)*weight_azimuthal(az_count,1)*fin_d(az_count,1)*((omega_y)^2*t_mki^3/3+omega_y*(y_old-Y_i_c(in_dx,in_dy))*t_mki^2+t_mki*(y_old-Y_i_c(in_dx,in_dy))^2);

                    M_xy_i(in_dx,in_dy)=M_xy_i(in_dx,in_dy)+1/(dx*dy)*mu(pol_count,1)*w(pol_count,1)*weight_azimuthal(az_count,1)*fin_d(az_count,1)*(omega_x*omega_y*t_mki^3/3+0.5*((x_old-X_i_c(in_dx,in_dy)*omega_y)+(y_old-Y_i_c(in_dx,in_dy))*omega_x)*t_mki^2+t_mki*(x_old-X_i_c(in_dx,in_dy))*(y_old-Y_i_c(in_dx,in_dy)));


                    x_new=dx*(in_dx-1);
                    y_new=tan(alt_azim_theta(az_count,1))*(x_new-x_old)+y_old;
        
                    if y_new<dy*(in_dy-1)&& abs(y_new-dy*(in_dy-1))>10^(-14)
                        y_new=dy*(in_dy-1);
                        x_new=x_old+(y_new-y_old)/tan(alt_azim_theta(az_count,1));
                         
                        in_dy=in_dy-1;
                        x_old=x_new;
                        y_old=y_new;
                     elseif abs(y_new-dy*(in_dy-1))<=10^(-14)
                        
                        in_dx=in_dx-1;
                        in_dy=in_dy-1;
                        x_old=x_new;
                        y_old=y_new;
                    else
                        
                        in_dx=in_dx-1;
                        x_old=x_new;
                        y_old=y_new;
                    end
                end
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
                    ray_index_count(in_dx,in_dy,az_count,pol_count)=ray_index_count(in_dx,in_dy,az_count,pol_count)+1;
                    omega_x=cos(alt_azim_theta(az_count,1))*mu(pol_count,1);
                    omega_y=sin(alt_azim_theta(az_count,1))*mu(pol_count,1);
                    x_c_track=x_c_t(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count,pol_count));
                    y_c_track=y_c_t(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count,pol_count));
                    t_mki=s_len(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))/ksi(in_dx,in_dy,az_count,1);
                    M_xx_i(in_dx,in_dy)=M_xx_i(in_dx,in_dy)+1/(dx*dy)*mu(pol_count,1)*w(pol_count,1)*weight_azimuthal(az_count,1)*fin_d(az_count,1)*((omega_x)^2*t_mki^3/3+omega_x*(x_old-X_i_c(in_dx,in_dy))*t_mki^2+t_mki*(x_old-X_i_c(in_dx,in_dy))^2);
                    M_yy_i(in_dx,in_dy)=M_yy_i(in_dx,in_dy)+1/(dx*dy)*mu(pol_count,1)*w(pol_count,1)*weight_azimuthal(az_count,1)*fin_d(az_count,1)*((omega_y)^2*t_mki^3/3+omega_y*(y_old-Y_i_c(in_dx,in_dy))*t_mki^2+t_mki*(y_old-Y_i_c(in_dx,in_dy))^2);

                    M_xy_i(in_dx,in_dy)=M_xy_i(in_dx,in_dy)+1/(dx*dy)*mu(pol_count,1)*w(pol_count,1)*weight_azimuthal(az_count,1)*fin_d(az_count,1)*(omega_x*omega_y*t_mki^3/3+0.5*((x_old-X_i_c(in_dx,in_dy)*omega_y)+(y_old-Y_i_c(in_dx,in_dy))*omega_x)*t_mki^2+t_mki*(x_old-X_i_c(in_dx,in_dy))*(y_old-Y_i_c(in_dx,in_dy)));

                    

                    x_new=dx*(in_dx-1);
                  
                    y_new=tan(alt_azim_theta(az_count,1))*(x_new-x_old)+y_old;
        
                    if y_new<dy*(in_dy-1) && abs(y_new-dy*(in_dy-1))>10^(-14)
                        y_new=dy*(in_dy-1);
                        x_new=x_old+(y_new-y_old)/tan(alt_azim_theta(az_count,1));
                        
                        in_dy=in_dy-1;
                        x_old=x_new;
                        y_old=y_new;
                     elseif abs(y_new-dy*(in_dy-1))<=10^(-14)
                        
                        in_dx=in_dx-1;
                        in_dy=in_dy-1;
                        x_old=x_new;
                        y_old=y_new;

                    else
                        
                        in_dx=in_dx-1;
                        x_old=x_new;
                        y_old=y_new;
                    end
                end
            end
    end

end

for i=1:mesh_center_ordinate_number
    for j=1:mesh_center_abscissa_number
        M(1,1,j,i)=1;
        M(2,2,j,i)=M_xx_i(j,i);
        M(2,3,j,i)=M_xy_i(j,i);
        M(3,2,j,i)=M_xy_i(j,i);
        M(3,3,j,i)=M_yy_i(j,i);
    end
end