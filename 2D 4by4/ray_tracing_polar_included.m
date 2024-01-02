function [exponential_portion,s_len,sum_s_len,altered_azimuthal_direction_theta,fin_d]=ray_tracing_polar_included

[altered_azimuthal_direction_theta,length_of_ray_in_a_mesh_to_particular_direction,fin_d]=ray_tracing;
%given data 

sigma_t=1;
sigma_s=0.7;
nu_sigma_f=0.39;


%spatial discretization
X=4;
Y=4;

dx=0.1;
dy=0.1;

x=(0:dx:X)';
y=(0:dy:Y)';

mesh_center_x=(dx/2:dx:X)';
mesh_center_y=(dy/2:dy:Y)';
mesh_center_abscissa_number=length(mesh_center_x);
mesh_center_ordinate_number=length(mesh_center_y);



%% angular discretization

%polar discretization

mu=[0.932954;0.537707;0.166648;-0.166648;-0.537707;-0.932954];
w=[0.670148;0.283619;0.046233;0.046233;0.283619;0.670148];
polar_discretization_number=size(mu,1);

%azimuthal discretization
N_a=128;
del_theta=2*pi/N_a;
theta=(0:del_theta:2*pi)';
azimuthal_direction_theta= 0.5*(theta(1:end-1,1)+theta(2:end,1));
azimuthal_discretization_number=size(azimuthal_direction_theta,1);

weight_azimuthal=zeros(N_a,1);
weight_azimuthal(2:end-1,1)=0.5*(altered_azimuthal_direction_theta(3:end,1)-altered_azimuthal_direction_theta(1:end-2,1));
weight_azimuthal(1,1)=0.5*(altered_azimuthal_direction_theta(1,1)+altered_azimuthal_direction_theta(2,1));
weight_azimuthal(end,1)=0.5*(2*pi-altered_azimuthal_direction_theta(end-1,1)+altered_azimuthal_direction_theta(1,1));

ray_index_count_for_each_mesh_for_each_direction=zeros(mesh_center_ordinate_number,mesh_center_abscissa_number,azimuthal_discretization_number,polar_discretization_number);
s_len=zeros(mesh_center_ordinate_number,mesh_center_abscissa_number,azimuthal_discretization_number,polar_discretization_number,100);
exponential_portion=zeros(mesh_center_ordinate_number,mesh_center_abscissa_number,azimuthal_discretization_number,polar_discretization_number,100);
sum_s_len=zeros(mesh_center_ordinate_number,mesh_center_abscissa_number,azimuthal_discretization_number,polar_discretization_number);

%% bottom to top rays

%left to right angle less than pi/2
for az_count=1:N_a/4
    %n_x=ceil(abs((X-dx/2)*sin(azimuthal_direction_theta(az_count,1))/init_d(az_count)));
    %n_y=ceil(abs((Y-dy/2)*cos(azimuthal_direction_theta(az_count,1))/init_d(az_count)));
    %altered_azimuthal_direction_theta(az_count,1)=atan(Y*n_x/(X*n_y));
    
    ray_spacing_x=fin_d(az_count,1)*abs(csc(altered_azimuthal_direction_theta(az_count,1)));
    ray_spacing_y=fin_d(az_count,1)*abs(sec(altered_azimuthal_direction_theta(az_count,1)));
    ray_pos_x_bound=(ray_spacing_x/2:ray_spacing_x:X)';
    ray_pos_y_bound=(ray_spacing_y/2:ray_spacing_y:Y)';


    number_of_ray_pos_x=length(ray_pos_x_bound);
    number_of_ray_pos_y=length(ray_pos_y_bound);
    for pol_count=1:polar_discretization_number
        p_x=1; % x index of ray position
    
      
        
        i_x=1; %mesh index x direction
        i_y=1; %mesh index y direction
    
        for p_x=1:number_of_ray_pos_x
            
            if(ray_pos_x_bound(p_x,1)>=dx*i_x)
                i_x=ceil(ray_pos_x_bound(p_x,1)/dx);
            end
            ray_tracking_index_x=i_x;
            ray_tracking_index_y=i_y;
            
            x_old=ray_pos_x_bound(p_x,1);
            y_old=0;
            
            while x_old<X && y_old<Y 
                x=dx*ray_tracking_index_x;
                y=tan(altered_azimuthal_direction_theta(az_count,1))*(x-x_old)+y_old;
                if y>dy*ray_tracking_index_y
                    ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count)=ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count)+1;
    
                    y_new=dy*ray_tracking_index_y;
                    x_new=(y_new-y_old)/tan(altered_azimuthal_direction_theta(az_count,1))+x_old;

                    s_len(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count, ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count))=length_of_ray_in_a_mesh_to_particular_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count))/abs(sin(mu(pol_count,1)));
                    exponential_portion(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count, ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count))=1-exp(-sigma_t*s_len(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count, ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count)));
                    sum_s_len(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count)=sum_s_len(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count)+ s_len(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count, ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count));
                    ray_tracking_index_y=ray_tracking_index_y+1;
                    
                    x_old=x_new;
                    y_old=y_new;
                else
                    ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count)=ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count)+1;
    
                    y_new=y;
                    x_new=x;
                    s_len(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count, ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count))=length_of_ray_in_a_mesh_to_particular_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count))/abs(sin(mu(pol_count,1)));


                    exponential_portion(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count, ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count))=1-exp(-sigma_t*s_len(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count, ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count)));
                    sum_s_len(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count)=sum_s_len(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count)+ s_len(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count, ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count));

                    x_old=x_new;
                    y_old=y_new;
                    ray_tracking_index_x=ray_tracking_index_x+1;
                    
                end
            end
        end
    
    
      
    
        p_y=1; % y index of ray position
        
        i_x=1; %mesh index x direction
        i_y=1; %mesh index y direction
        for p_y=1:number_of_ray_pos_y
            
            if(ray_pos_y_bound(p_y,1)>=dy*i_y)
                i_y=ceil(ray_pos_y_bound(p_y,1)/dy);
            end
            ray_tracking_index_x=i_x;
            ray_tracking_index_y=i_y;
            
            x_old=0;
            y_old=ray_pos_y_bound(p_y,1);
            
            while x_old<X && y_old<Y 
                x=dx*ray_tracking_index_x;
                y=tan(altered_azimuthal_direction_theta(az_count,1))*(x-x_old)+y_old;
                if y>dy*ray_tracking_index_y
                    ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count)=ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count)+1;
    
                    y_new=dy*ray_tracking_index_y;
                    x_new=(y_new-y_old)/tan(altered_azimuthal_direction_theta(az_count,1))+x_old;
                    s_len(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count, ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count))=length_of_ray_in_a_mesh_to_particular_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count))/abs(sin(mu(pol_count,1)));


                    exponential_portion(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count, ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count))=1-exp(-sigma_t*s_len(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count, ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count)));
                    sum_s_len(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count)=sum_s_len(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count)+ s_len(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count, ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count));

                    ray_tracking_index_y=ray_tracking_index_y+1;
                   
                    x_old=x_new;
                    y_old=y_new;
                else
                    ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count)=ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count)+1;

                    y_new=y;
                    x_new=x;
                    s_len(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count, ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count))=length_of_ray_in_a_mesh_to_particular_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count))/abs(sin(mu(pol_count,1)));


                    exponential_portion(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count, ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count))=1-exp(-sigma_t*s_len(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count, ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count)));
                    sum_s_len(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count)=sum_s_len(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count)+ s_len(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count, ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count));

                    x_old=x_new;
                    y_old=y_new;
                    ray_tracking_index_x=ray_tracking_index_x+1;
                    
                end
            end
        end
    end
end


%right to left angle greater than pi/2
for az_count=N_a/4+1:N_a/2
    %n_x=ceil(abs((X-dx/2)*sin(azimuthal_direction_theta(az_count,1))/init_d(az_count)));
    %n_y=ceil(abs((Y-dy/2)*cos(azimuthal_direction_theta(az_count,1))/init_d(az_count)));
    %altered_azimuthal_direction_theta(az_count,1)=pi-atan(Y*n_x/(X*n_y));
    %fin_d(az_count,1)=sqrt((X^2+Y^2)/(n_x^2+n_y^2));

    ray_spacing_x=fin_d(az_count,1)*abs(csc(altered_azimuthal_direction_theta(az_count,1)));
    ray_spacing_y=fin_d(az_count,1)*abs(sec(altered_azimuthal_direction_theta(az_count,1)));
    ray_pos_x_bound=(ray_spacing_x/2:ray_spacing_x:X)';
    ray_pos_y_bound=(ray_spacing_y/2:ray_spacing_y:Y)';


    number_of_ray_pos_x=length(ray_pos_x_bound);
    number_of_ray_pos_y=length(ray_pos_y_bound);
    for pol_count=1:polar_discretization_number

        p_x=number_of_ray_pos_x; % x index of ray position
    
        p_y=1; % y index of ray position
        
        i_x=mesh_center_abscissa_number; %mesh index x direction
        i_y=1; %mesh index y direction
    
        for p_x=number_of_ray_pos_x:-1:1
            
            %if(ray_pos_x_bound(p_x,1)<=dx*i_x)
                i_x=ceil(ray_pos_x_bound(p_x,1)/dx);
            %end
            ray_tracking_index_x=i_x;
            ray_tracking_index_y=i_y;
            
            x_old=ray_pos_x_bound(p_x,1);
            y_old=0;
            
            while x_old>0 && y_old<Y 
                x=dx*(ray_tracking_index_x-1);
                y=tan(altered_azimuthal_direction_theta(az_count,1))*(x-x_old)+y_old;
                if y>dy*ray_tracking_index_y
                    ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count)=ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count)+1;

                    y_new=dy*ray_tracking_index_y;
                    x_new=(y_new-y_old)/tan(altered_azimuthal_direction_theta(az_count,1))+x_old;
                    s_len(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count, ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count))=length_of_ray_in_a_mesh_to_particular_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count))/abs(sin(mu(pol_count,1)));


                    exponential_portion(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count, ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count))=1-exp(-sigma_t*s_len(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count, ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count)));
                    sum_s_len(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count)=sum_s_len(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count)+ s_len(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count, ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count));

                    ray_tracking_index_y=ray_tracking_index_y+1;
                    
                    x_old=x_new;
                    y_old=y_new;
                else
                    ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count)=ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count)+1;
                    
                    y_new=y;
                    x_new=x;
                   s_len(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count, ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count))=length_of_ray_in_a_mesh_to_particular_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count))/abs(sin(mu(pol_count,1)));


                    exponential_portion(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count, ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count))=1-exp(-sigma_t*s_len(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count, ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count)));
                    sum_s_len(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count)=sum_s_len(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count)+ s_len(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count, ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count));

                    x_old=x_new;
                    y_old=y_new;
                    ray_tracking_index_x=ray_tracking_index_x-1;
                    
                end
            end
        end
        
        p_y=1; % y index of ray position
        
        i_x=mesh_center_abscissa_number; %mesh index x direction
        i_y=1; %mesh index y direction
    
        for p_y=1:number_of_ray_pos_y
            
            if(ray_pos_y_bound(p_y,1)>=dy*i_y)
                i_y=ceil(ray_pos_y_bound(p_y,1)/dy);
            end
            ray_tracking_index_x=i_x;
            ray_tracking_index_y=i_y;
            
            x_old=X;
            y_old=ray_pos_y_bound(p_y,1);
            
            while x_old>0 && y_old<Y 
                x=dx*(ray_tracking_index_x-1);
                y=tan(altered_azimuthal_direction_theta(az_count,1))*(x-x_old)+y_old;
                if y>dy*ray_tracking_index_y
                    ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count)=ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count)+1;

                    y_new=dy*ray_tracking_index_y;
                    x_new=(y_new-y_old)/tan(altered_azimuthal_direction_theta(az_count,1))+x_old;
                    s_len(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count, ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count))=length_of_ray_in_a_mesh_to_particular_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count))/abs(sin(mu(pol_count,1)));


                    exponential_portion(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count, ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count))=1-exp(-sigma_t*s_len(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count, ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count)));
                    sum_s_len(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count)=sum_s_len(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count)+ s_len(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count, ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count));

                    ray_tracking_index_y=ray_tracking_index_y+1;
                    
                    x_old=x_new;
                    y_old=y_new;
                else
                    ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count)=ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count)+1;

                    y_new=y;
                    x_new=x;
                   s_len(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count, ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count))=length_of_ray_in_a_mesh_to_particular_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count))/abs(sin(mu(pol_count,1)));


                    exponential_portion(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count, ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count))=1-exp(-sigma_t*s_len(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count, ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count)));
                    sum_s_len(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count)=sum_s_len(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count)+ s_len(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count, ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count));

                    x_old=x_new;
                    y_old=y_new;
                    ray_tracking_index_x=ray_tracking_index_x-1;
                    
                end
            end
        end
    end
end



%% top to bottom rays

%left to right angle less than 2*pi greater than 3*pi/2
for az_count=3*N_a/4+1:N_a
    %n_x=ceil(abs((X-dx/2)*sin(azimuthal_direction_theta(az_count,1))/init_d(az_count)));
   % n_y=ceil(abs((Y-dy/2)*cos(azimuthal_direction_theta(az_count,1))/init_d(az_count)));
    %altered_azimuthal_direction_theta(az_count,1)=2*pi-atan(Y*n_x/(X*n_y));
  

    ray_spacing_x=fin_d(az_count,1)*abs(csc(altered_azimuthal_direction_theta(az_count,1)));
    ray_spacing_y=fin_d(az_count,1)*abs(sec(altered_azimuthal_direction_theta(az_count,1)));
    ray_pos_x_bound=(ray_spacing_x/2:ray_spacing_x:X)';
    ray_pos_y_bound=(ray_spacing_y/2:ray_spacing_y:Y)';

    number_of_ray_pos_x=length(ray_pos_x_bound);
    number_of_ray_pos_y=length(ray_pos_y_bound);
    for pol_count=1:polar_discretization_number

        p_x=1; % x index of ray position
    
      
       
        i_x=1; %mesh index x direction
        i_y=mesh_center_ordinate_number; %mesh index y direction
    
        for p_x=1:number_of_ray_pos_x
            
            if(ray_pos_x_bound(p_x,1)>=dx*i_x)
                i_x=ceil(ray_pos_x_bound(p_x,1)/dx);
            end
            ray_tracking_index_x=i_x;
            ray_tracking_index_y=i_y;
            
            x_old=ray_pos_x_bound(p_x,1);
            y_old=Y;
            
            while x_old<X && y_old>0 
                x=dx*ray_tracking_index_x;
                y=tan(altered_azimuthal_direction_theta(az_count,1))*(x-x_old)+y_old;
                if y<dy*(ray_tracking_index_y-1)
                    ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count)=ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count)+1;

                    y_new=dy*(ray_tracking_index_y-1);
                    x_new=(y_new-y_old)/tan(altered_azimuthal_direction_theta(az_count,1))+x_old;
                  s_len(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count, ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count))=length_of_ray_in_a_mesh_to_particular_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count))/abs(sin(mu(pol_count,1)));


                    exponential_portion(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count, ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count))=1-exp(-sigma_t*s_len(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count, ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count)));
                    sum_s_len(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count)=sum_s_len(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count)+ s_len(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count, ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count));

                    ray_tracking_index_y=ray_tracking_index_y-1;
                    
                    x_old=x_new;
                    y_old=y_new;
                else
                    ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count)=ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count)+1;

                    y_new=y;
                    x_new=x;
                  s_len(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count, ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count))=length_of_ray_in_a_mesh_to_particular_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count))/abs(sin(mu(pol_count,1)));


                    exponential_portion(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count, ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count))=1-exp(-sigma_t*s_len(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count, ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count)));
                    sum_s_len(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count)=sum_s_len(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count)+ s_len(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count, ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count));

                    x_old=x_new;
                    y_old=y_new;
                    ray_tracking_index_x=ray_tracking_index_x+1;
                    
                end
            end
        end
    
    
    
    
        p_y=1; % y index of ray position
        
        i_x=1; %mesh index x direction
        i_y=mesh_center_ordinate_number; %mesh index y direction
        for p_y=number_of_ray_pos_y:-1:1
            
            %if(ray_pos_y_bound(p_y,1)>=dy*i_y)
                i_y=ceil(ray_pos_y_bound(p_y,1)/dy);
            %end
            ray_tracking_index_x=i_x;
            ray_tracking_index_y=i_y;
            
            x_old=0;
            y_old=ray_pos_y_bound(p_y,1);
            
            while x_old<X && y_old>0 
                x=dx*ray_tracking_index_x;
                y=tan(altered_azimuthal_direction_theta(az_count,1))*(x-x_old)+y_old;
                if y<dy*(ray_tracking_index_y-1)
                    ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count)=ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count)+1;

                    y_new=dy*(ray_tracking_index_y-1);
                    x_new=(y_new-y_old)/tan(altered_azimuthal_direction_theta(az_count,1))+x_old;
                    s_len(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count, ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count))=length_of_ray_in_a_mesh_to_particular_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count))/abs(sin(mu(pol_count,1)));


                    exponential_portion(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count, ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count))=1-exp(-sigma_t*s_len(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count, ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count)));
                    sum_s_len(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count)=sum_s_len(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count)+ s_len(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count, ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count));

                    ray_tracking_index_y=ray_tracking_index_y-1;
                   
                    x_old=x_new;
                    y_old=y_new;
                else
                    ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count)=ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count)+1;

                    y_new=y;
                    x_new=x;
                   s_len(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count, ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count))=length_of_ray_in_a_mesh_to_particular_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count))/abs(sin(mu(pol_count,1)));


                    exponential_portion(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count, ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count))=1-exp(-sigma_t*s_len(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count, ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count)));
                    sum_s_len(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count)=sum_s_len(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count)+ s_len(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count, ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count));

                    x_old=x_new;
                    y_old=y_new;
                    ray_tracking_index_x=ray_tracking_index_x+1;
                    
                end
            end
        end
    end

end
  

%right to left angle greater than pi/2
for az_count=N_a/2+1:3*N_a/4
    %n_x=ceil(abs((X-dx/2)*sin(azimuthal_direction_theta(az_count,1))/init_d(az_count)));
    %n_y=ceil(abs((Y-dy/2)*cos(azimuthal_direction_theta(az_count,1))/init_d(az_count)));
    %altered_azimuthal_direction_theta(az_count,1)=pi+atan(Y*n_x/(X*n_y));
 
    ray_spacing_x=fin_d(az_count,1)*abs(csc(altered_azimuthal_direction_theta(az_count,1)));
    ray_spacing_y=fin_d(az_count,1)*abs(sec(altered_azimuthal_direction_theta(az_count,1)));
   ray_pos_x_bound=(ray_spacing_x/2:ray_spacing_x:X)';
    ray_pos_y_bound=(ray_spacing_y/2:ray_spacing_y:Y)';


    number_of_ray_pos_x=length(ray_pos_x_bound);
    number_of_ray_pos_y=length(ray_pos_y_bound);
    for pol_count=1:polar_discretization_number

        p_x=number_of_ray_pos_x; % x index of ray position
    
        p_y=number_of_ray_pos_y; % y index of ray position
        
        i_x=mesh_center_abscissa_number; %mesh index x direction
        i_y=mesh_center_ordinate_number; %mesh index y direction
    
        for p_x=number_of_ray_pos_x:-1:1
            
            %if(ray_pos_x_bound(p_x,1)<=dx*i_x)
                i_x=ceil(ray_pos_x_bound(p_x,1)/dx);
            %end
            ray_tracking_index_x=i_x;
            ray_tracking_index_y=i_y;
            
            x_old=ray_pos_x_bound(p_x,1);
            y_old=Y;
            
            while x_old>0 && y_old>0
                x=dx*(ray_tracking_index_x-1);
                y=tan(altered_azimuthal_direction_theta(az_count,1))*(x-x_old)+y_old;
                if y<dy*(ray_tracking_index_y-1)
                    ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count)=ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count)+1;

                    y_new=dy*(ray_tracking_index_y-1);
                    x_new=(y_new-y_old)/tan(altered_azimuthal_direction_theta(az_count,1))+x_old;
                    s_len(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count, ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count))=length_of_ray_in_a_mesh_to_particular_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count))/abs(sin(mu(pol_count,1)));


                    exponential_portion(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count, ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count))=1-exp(-sigma_t*s_len(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count, ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count)));
                    sum_s_len(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count)=sum_s_len(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count)+ s_len(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count, ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count));

                    ray_tracking_index_y=ray_tracking_index_y-1;
                   
                    x_old=x_new;
                    y_old=y_new;
                else
                    ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count)=ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count)+1;

                    y_new=y;
                    x_new=x;
                   s_len(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count, ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count))=length_of_ray_in_a_mesh_to_particular_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count))/abs(sin(mu(pol_count,1)));


                    exponential_portion(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count, ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count))=1-exp(-sigma_t*s_len(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count, ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count)));
                    sum_s_len(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count)=sum_s_len(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count)+ s_len(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count, ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count));

                    x_old=x_new;
                    y_old=y_new;
                    ray_tracking_index_x=ray_tracking_index_x-1;
                    
                end
            end
        end
        
        p_y=1; % y index of ray position
        
        i_x=mesh_center_abscissa_number; %mesh index x direction
        i_y=mesh_center_ordinate_number; %mesh index y direction
    
        for p_y=number_of_ray_pos_y:-1:1
            
            %if(ray_pos_y_bound(p_y,1)>=dy*i_y)
                i_y=ceil(ray_pos_y_bound(p_y,1)/dy);
            %end
            ray_tracking_index_x=i_x;
            ray_tracking_index_y=i_y;
            
            x_old=X;
            y_old=ray_pos_y_bound(p_y,1);
            
            while x_old>0 && y_old>0 
                x=dx*(ray_tracking_index_x-1);
                y=tan(altered_azimuthal_direction_theta(az_count,1))*(x-x_old)+y_old;
                if y<dy*(ray_tracking_index_y-1)
                    ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count)=ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count)+1;

                    y_new=dy*(ray_tracking_index_y-1);
                    x_new=(y_new-y_old)/tan(altered_azimuthal_direction_theta(az_count,1))+x_old;
                    s_len(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count, ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count))=length_of_ray_in_a_mesh_to_particular_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count))/abs(sin(mu(pol_count,1)));


                    exponential_portion(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count, ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count))=1-exp(-sigma_t*s_len(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count, ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count)));
                    sum_s_len(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count)=sum_s_len(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count)+ s_len(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count, ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count));

                    ray_tracking_index_y=ray_tracking_index_y-1;
                   
                    x_old=x_new;
                    y_old=y_new;
                else
                    ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count)=ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count)+1;
                    
                    y_new=y;
                    x_new=x;
                    s_len(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count, ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count))=length_of_ray_in_a_mesh_to_particular_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count))/abs(sin(mu(pol_count,1)));
                    

                    exponential_portion(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count, ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count))=1-exp(-sigma_t*s_len(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count, ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count)));
                    sum_s_len(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count)=sum_s_len(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count)+ s_len(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count, ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,pol_count));
                    
                    x_old=x_new;
                    y_old=y_new;
                    ray_tracking_index_x=ray_tracking_index_x-1;
                    
                end
            end
        end
    end
end


%% sweep code




