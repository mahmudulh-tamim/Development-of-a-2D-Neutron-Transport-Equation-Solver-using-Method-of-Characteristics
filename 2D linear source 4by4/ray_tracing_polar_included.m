
function [F_1,F_2,G_1,G_2,H,tau,ksi,x_c_t,y_c_t,X_i_c,Y_i_c,s_len,sum_s_len,adj_len,alt_azim_theta,fin_d]=ray_tracing_polar_included(X,Y,dx,dy,N_a,sigma_t)


[alt_azim_theta,length_of_rays,fin_d,sum_len,x_c_t,y_c_t]=ray_tracing(X,Y,dx,dy,N_a);





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


s_len=zeros(mesh_center_abscissa_number,mesh_center_ordinate_number,azimuthal_discretization_number,polar_discretization_number,500);

tau=zeros(mesh_center_abscissa_number,mesh_center_ordinate_number,azimuthal_discretization_number,polar_discretization_number,500);

F_1=zeros(mesh_center_abscissa_number,mesh_center_ordinate_number,azimuthal_discretization_number,polar_discretization_number,500);
F_2=zeros(mesh_center_abscissa_number,mesh_center_ordinate_number,azimuthal_discretization_number,polar_discretization_number,500);
G_1=zeros(mesh_center_abscissa_number,mesh_center_ordinate_number,azimuthal_discretization_number,polar_discretization_number,500);
G_2=zeros(mesh_center_abscissa_number,mesh_center_ordinate_number,azimuthal_discretization_number,polar_discretization_number,500);
H=zeros(mesh_center_abscissa_number,mesh_center_ordinate_number,azimuthal_discretization_number,polar_discretization_number,500);
X_i_c=zeros(mesh_center_abscissa_number,mesh_center_ordinate_number);
Y_i_c=zeros(mesh_center_abscissa_number,mesh_center_ordinate_number);


M=zeros(3,3,mesh_center_abscissa_number,mesh_center_ordinate_number);


sum_s_len=zeros(mesh_center_abscissa_number,mesh_center_ordinate_number,azimuthal_discretization_number,polar_discretization_number);
sum_s_len_red=zeros(mesh_center_abscissa_number,mesh_center_ordinate_number,azimuthal_discretization_number,polar_discretization_number);
adj_len=length_of_rays;
area_approx=zeros(mesh_center_abscissa_number,mesh_center_ordinate_number,azimuthal_discretization_number);

total_rays=zeros(azimuthal_discretization_number,polar_discretization_number);

for azim=1:N_a
    area_approx(:,:,azim)=sum_len(:,:,azim)*fin_d(azim,1);
end
ksi=dx*dy./area_approx;

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
                    adj_len(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count,pol_count))=length_of_rays(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count,pol_count))*ksi(in_dx,in_dy,az_count);

                    s_len(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))=adj_len(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count,pol_count))/abs(mu(pol_count,1));
                    tau(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))= s_len(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))*sigma_t;

                    F_1(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))=1-exp(-tau(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count)));
                    F_2(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))=2*(tau(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))-F_1(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count)))-tau(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))*F_1(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count));

                    G_1(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))=1+tau(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))/2-(1+1/tau(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count)))*F_1(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count));
                    G_2(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))=(2/3)*tau(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))-(1+2/tau(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count)))*G_1(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count));
                    H(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))=tau(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))/2-G_1(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count));
                    X_i_c(in_dx,in_dy)=X_i_c(in_dx,in_dy)+mu(pol_count,1)*s_len(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))*x_c_t(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count,pol_count))*w(pol_count,1)*weight_azimuthal(az_count,1)*fin_d(az_count,1)/(ksi(in_dx,in_dy,az_count)*dx*dy*4*pi);
                    Y_i_c(in_dx,in_dy)=Y_i_c(in_dx,in_dy)+mu(pol_count,1)*s_len(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))*y_c_t(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count,pol_count))*w(pol_count,1)*weight_azimuthal(az_count,1)*fin_d(az_count,1)/(ksi(in_dx,in_dy,az_count)*dx*dy*4*pi);
                    
                    
                    sum_s_len(in_dx,in_dy,az_count,pol_count)=sum_s_len(in_dx,in_dy,az_count,pol_count)+s_len(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count));
                    sum_s_len_red(in_dx,in_dy,az_count,pol_count)=sum_len(in_dx,in_dy,az_count)/abs(mu(pol_count,1));


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
                    adj_len(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count,pol_count))=length_of_rays(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count,pol_count))*ksi(in_dx,in_dy,az_count);

                    s_len(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))=adj_len(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count,pol_count))/abs(mu(pol_count,1));
                    tau(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))= s_len(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))*sigma_t;

                    F_1(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))=1-exp(-tau(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count)));
                    F_2(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))=2*(tau(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))-F_1(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count)))-tau(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))*F_1(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count));

                    G_1(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))=1+tau(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))/2-(1+1/tau(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count)))*F_1(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count));
                    G_2(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))=(2/3)*tau(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))-(1+2/tau(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count)))*G_1(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count));
                    H(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))=tau(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))/2-G_1(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count));

                    X_i_c(in_dx,in_dy)=X_i_c(in_dx,in_dy)+mu(pol_count,1)*s_len(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))*x_c_t(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count,pol_count))*w(pol_count,1)*weight_azimuthal(az_count,1)*fin_d(az_count,1)/(ksi(in_dx,in_dy,az_count)*dx*dy*4*pi);
                    Y_i_c(in_dx,in_dy)=Y_i_c(in_dx,in_dy)+mu(pol_count,1)*s_len(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))*y_c_t(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count,pol_count))*w(pol_count,1)*weight_azimuthal(az_count,1)*fin_d(az_count,1)/(ksi(in_dx,in_dy,az_count)*dx*dy*4*pi);
                    
                    
                    sum_s_len(in_dx,in_dy,az_count,pol_count)=sum_s_len(in_dx,in_dy,az_count,pol_count)+s_len(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count));
                    sum_s_len_red(in_dx,in_dy,az_count,pol_count)=sum_len(in_dx,in_dy,az_count)/abs(mu(pol_count,1));

                    
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
                    adj_len(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count,pol_count))=length_of_rays(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count,pol_count))*ksi(in_dx,in_dy,az_count);

                    s_len(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))=adj_len(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count,pol_count))/abs(mu(pol_count,1));
                    tau(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))= s_len(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))*sigma_t;

                    F_1(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))=1-exp(-tau(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count)));
                    F_2(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))=2*(tau(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))-F_1(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count)))-tau(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))*F_1(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count));

                    G_1(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))=1+tau(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))/2-(1+1/tau(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count)))*F_1(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count));
                    G_2(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))=(2/3)*tau(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))-(1+2/tau(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count)))*G_1(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count));
                    H(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))=tau(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))/2-G_1(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count));

                    X_i_c(in_dx,in_dy)=X_i_c(in_dx,in_dy)+mu(pol_count,1)*s_len(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))*x_c_t(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count,pol_count))*w(pol_count,1)*weight_azimuthal(az_count,1)*fin_d(az_count,1)/(ksi(in_dx,in_dy,az_count)*dx*dy*4*pi);
                    Y_i_c(in_dx,in_dy)=Y_i_c(in_dx,in_dy)+mu(pol_count,1)*s_len(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))*y_c_t(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count,pol_count))*w(pol_count,1)*weight_azimuthal(az_count,1)*fin_d(az_count,1)/(ksi(in_dx,in_dy,az_count)*dx*dy*4*pi);
                    
                    
                    sum_s_len(in_dx,in_dy,az_count,pol_count)=sum_s_len(in_dx,in_dy,az_count,pol_count)+s_len(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count));
                    sum_s_len_red(in_dx,in_dy,az_count,pol_count)=sum_len(in_dx,in_dy,az_count)/abs(mu(pol_count,1));


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
                    adj_len(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count,pol_count))=length_of_rays(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count,pol_count))*ksi(in_dx,in_dy,az_count);

                    s_len(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))=adj_len(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count,pol_count))/abs(mu(pol_count,1));
                    tau(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))= s_len(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))*sigma_t;

                    F_1(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))=1-exp(-tau(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count)));
                    F_2(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))=2*(tau(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))-F_1(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count)))-tau(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))*F_1(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count));

                    G_1(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))=1+tau(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))/2-(1+1/tau(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count)))*F_1(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count));
                    G_2(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))=(2/3)*tau(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))-(1+2/tau(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count)))*G_1(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count));
                    H(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))=tau(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))/2-G_1(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count));

                    X_i_c(in_dx,in_dy)=X_i_c(in_dx,in_dy)+mu(pol_count,1)*s_len(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))*x_c_t(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count,pol_count))*w(pol_count,1)*weight_azimuthal(az_count,1)*fin_d(az_count,1)/(ksi(in_dx,in_dy,az_count)*dx*dy*4*pi);
                    Y_i_c(in_dx,in_dy)=Y_i_c(in_dx,in_dy)+mu(pol_count,1)*s_len(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))*y_c_t(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count,pol_count))*w(pol_count,1)*weight_azimuthal(az_count,1)*fin_d(az_count,1)/(ksi(in_dx,in_dy,az_count)*dx*dy*4*pi);
                    
                    
                    sum_s_len(in_dx,in_dy,az_count,pol_count)=sum_s_len(in_dx,in_dy,az_count,pol_count)+s_len(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count));
                    sum_s_len_red(in_dx,in_dy,az_count,pol_count)=sum_len(in_dx,in_dy,az_count)/abs(mu(pol_count,1));


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
                    
                    adj_len(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count,pol_count))=length_of_rays(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count,pol_count))*ksi(in_dx,in_dy,az_count);

                    s_len(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))=adj_len(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count,pol_count))/abs(mu(pol_count,1));
                    tau(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))= s_len(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))*sigma_t;

                    F_1(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))=1-exp(-tau(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count)));
                    F_2(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))=2*(tau(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))-F_1(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count)))-tau(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))*F_1(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count));

                    G_1(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))=1+tau(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))/2-(1+1/tau(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count)))*F_1(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count));
                    G_2(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))=(2/3)*tau(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))-(1+2/tau(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count)))*G_1(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count));
                    H(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))=tau(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))/2-G_1(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count));

                    X_i_c(in_dx,in_dy)=X_i_c(in_dx,in_dy)+mu(pol_count,1)*s_len(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))*x_c_t(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count,pol_count))*w(pol_count,1)*weight_azimuthal(az_count,1)*fin_d(az_count,1)/(ksi(in_dx,in_dy,az_count)*dx*dy*4*pi);
                    Y_i_c(in_dx,in_dy)=Y_i_c(in_dx,in_dy)+mu(pol_count,1)*s_len(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))*y_c_t(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count,pol_count))*w(pol_count,1)*weight_azimuthal(az_count,1)*fin_d(az_count,1)/(ksi(in_dx,in_dy,az_count)*dx*dy*4*pi);
                    
                    
                    sum_s_len(in_dx,in_dy,az_count,pol_count)=sum_s_len(in_dx,in_dy,az_count,pol_count)+s_len(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count));
                    sum_s_len_red(in_dx,in_dy,az_count,pol_count)=sum_len(in_dx,in_dy,az_count)/abs(mu(pol_count,1));


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

                    adj_len(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count,pol_count))=length_of_rays(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count,pol_count))*ksi(in_dx,in_dy,az_count);

                    s_len(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))=adj_len(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count,pol_count))/abs(mu(pol_count,1));
                    tau(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))= s_len(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))*sigma_t;

                    F_1(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))=1-exp(-tau(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count)));
                    F_2(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))=2*(tau(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))-F_1(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count)))-tau(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))*F_1(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count));

                    G_1(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))=1+tau(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))/2-(1+1/tau(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count)))*F_1(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count));
                    G_2(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))=(2/3)*tau(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))-(1+2/tau(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count)))*G_1(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count));
                    H(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))=tau(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))/2-G_1(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count));

                    X_i_c(in_dx,in_dy)=X_i_c(in_dx,in_dy)+mu(pol_count,1)*s_len(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))*x_c_t(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count,pol_count))*w(pol_count,1)*weight_azimuthal(az_count,1)*fin_d(az_count,1)/(ksi(in_dx,in_dy,az_count)*dx*dy*4*pi);
                    Y_i_c(in_dx,in_dy)=Y_i_c(in_dx,in_dy)+mu(pol_count,1)*s_len(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))*y_c_t(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count,pol_count))*w(pol_count,1)*weight_azimuthal(az_count,1)*fin_d(az_count,1)/(ksi(in_dx,in_dy,az_count)*dx*dy*4*pi);
                    
                    
                    sum_s_len(in_dx,in_dy,az_count,pol_count)=sum_s_len(in_dx,in_dy,az_count,pol_count)+s_len(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count));
                    sum_s_len_red(in_dx,in_dy,az_count,pol_count)=sum_len(in_dx,in_dy,az_count)/abs(mu(pol_count,1));


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
                    adj_len(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count,pol_count))=length_of_rays(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count,pol_count))*ksi(in_dx,in_dy,az_count);

                    s_len(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))=adj_len(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count,pol_count))/abs(mu(pol_count,1));
                    tau(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))= s_len(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))*sigma_t;

                    F_1(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))=1-exp(-tau(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count)));
                    F_2(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))=2*(tau(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))-F_1(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count)))-tau(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))*F_1(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count));

                    G_1(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))=1+tau(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))/2-(1+1/tau(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count)))*F_1(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count));
                    G_2(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))=(2/3)*tau(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))-(1+2/tau(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count)))*G_1(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count));
                    H(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))=tau(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))/2-G_1(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count));

                    X_i_c(in_dx,in_dy)=X_i_c(in_dx,in_dy)+mu(pol_count,1)*s_len(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))*x_c_t(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count,pol_count))*w(pol_count,1)*weight_azimuthal(az_count,1)*fin_d(az_count,1)/(ksi(in_dx,in_dy,az_count)*dx*dy*4*pi);
                    Y_i_c(in_dx,in_dy)=Y_i_c(in_dx,in_dy)+mu(pol_count,1)*s_len(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))*y_c_t(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count,pol_count))*w(pol_count,1)*weight_azimuthal(az_count,1)*fin_d(az_count,1)/(ksi(in_dx,in_dy,az_count)*dx*dy*4*pi);
                    
                    
                    sum_s_len(in_dx,in_dy,az_count,pol_count)=sum_s_len(in_dx,in_dy,az_count,pol_count)+s_len(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count));
                    sum_s_len_red(in_dx,in_dy,az_count,pol_count)=sum_len(in_dx,in_dy,az_count)/abs(mu(pol_count,1));


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
                    
                    adj_len(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count,pol_count))=length_of_rays(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count,pol_count))*ksi(in_dx,in_dy,az_count);

                    s_len(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))=adj_len(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count,pol_count))/abs(mu(pol_count,1));
                    tau(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))= s_len(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))*sigma_t;

                    F_1(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))=1-exp(-tau(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count)));
                    F_2(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))=2*(tau(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))-F_1(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count)))-tau(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))*F_1(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count));

                    G_1(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))=1+tau(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))/2-(1+1/tau(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count)))*F_1(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count));
                    G_2(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))=(2/3)*tau(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))-(1+2/tau(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count)))*G_1(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count));
                    H(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))=tau(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))/2-G_1(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count));
                    
                    X_i_c(in_dx,in_dy)=X_i_c(in_dx,in_dy)+mu(pol_count,1)*s_len(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))*x_c_t(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count,pol_count))*w(pol_count,1)*weight_azimuthal(az_count,1)*fin_d(az_count,1)/(ksi(in_dx,in_dy,az_count)*dx*dy*4*pi);
                    Y_i_c(in_dx,in_dy)=Y_i_c(in_dx,in_dy)+mu(pol_count,1)*s_len(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count))*y_c_t(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count,pol_count))*w(pol_count,1)*weight_azimuthal(az_count,1)*fin_d(az_count,1)/(ksi(in_dx,in_dy,az_count)*dx*dy*4*pi);
                    
                    
                    sum_s_len(in_dx,in_dy,az_count,pol_count)=sum_s_len(in_dx,in_dy,az_count,pol_count)+s_len(in_dx,in_dy,az_count,pol_count,ray_index_count(in_dx,in_dy,az_count,pol_count));
                    sum_s_len_red(in_dx,in_dy,az_count,pol_count)=sum_len(in_dx,in_dy,az_count)/abs(mu(pol_count,1));


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

