%% bottom to top rays

%left to right angle less than pi/2
for az_count=1:N_a/4
    n_x=ceil(abs((X-dx/2)*sin(azimuthal_direction_theta(az_count,1))/init_d(az_count)));
    n_y=ceil(abs((Y-dy/2)*cos(azimuthal_direction_theta(az_count,1))/init_d(az_count)));
    altered_azimuthal_direction_theta(az_count,1)=atan(Y*n_x/(X*n_y));
    fin_d(az_count,1)=sqrt((X^2+Y^2)/(n_x^2+n_y^2));
    ray_spacing_x=fin_d(az_count,1)*abs(csc(altered_azimuthal_direction_theta(az_count,1)));
    ray_spacing_y=fin_d(az_count,1)*abs(sec(altered_azimuthal_direction_theta(az_count,1)));
    ray_pos_x_bound=(dx/4:ray_spacing_x:X-dx/4)';
    ray_pos_y_bound=(dy/4:ray_spacing_y:Y-dy/4)';


    number_of_ray_pos_x=length(ray_pos_x_bound);
    number_of_ray_pos_y=length(ray_pos_y_bound);
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
                ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count)=ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count)+1;

                y_new=dy*ray_tracking_index_y;
                x_new=(y_new-y_old)/tan(altered_azimuthal_direction_theta(az_count,1))+x_old;
                length_of_ray_in_a_mesh_to_particular_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count))=sqrt((x_new-x_old)^2+(y_new-y_old)^2);
                ray_tracking_index_y=ray_tracking_index_y+1;
                plot([x_old, x_new], [y_old, y_new]);
                x_old=x_new;
                y_old=y_new;
            else
                ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count)=ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count)+1;

                y_new=y;
                x_new=x;
                length_of_ray_in_a_mesh_to_particular_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count))=sqrt((x_new-x_old)^2+(y_new-y_old)^2);
                plot([x_old, x_new], [y_old, y_new]);
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
                ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count)=ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count)+1;

                y_new=dy*ray_tracking_index_y;
                x_new=(y_new-y_old)/tan(altered_azimuthal_direction_theta(az_count,1))+x_old;
                length_of_ray_in_a_mesh_to_particular_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count))=sqrt((x_new-x_old)^2+(y_new-y_old)^2);
                ray_tracking_index_y=ray_tracking_index_y+1;
                plot([x_old, x_new], [y_old, y_new]);
                x_old=x_new;
                y_old=y_new;
            else
                ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count)=ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count)+1;

                y_new=y;
                x_new=x;
                length_of_ray_in_a_mesh_to_particular_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count))=sqrt((x_new-x_old)^2+(y_new-y_old)^2);
                plot([x_old, x_new], [y_old, y_new]);
                x_old=x_new;
                y_old=y_new;
                ray_tracking_index_x=ray_tracking_index_x+1;
                
            end
        end
    end
end


%right to left angle greater than pi/2
for az_count=N_a/4+1:N_a/2
    n_x=ceil(abs((X-dx/2)*sin(azimuthal_direction_theta(az_count,1))/init_d(az_count)));
    n_y=ceil(abs((Y-dy/2)*cos(azimuthal_direction_theta(az_count,1))/init_d(az_count)));
    altered_azimuthal_direction_theta(az_count,1)=pi-atan(Y*n_x/(X*n_y));
    fin_d(az_count,1)=sqrt((X^2+Y^2)/(n_x^2+n_y^2));
    ray_spacing_x=fin_d(az_count,1)*abs(csc(altered_azimuthal_direction_theta(az_count,1)));
    ray_spacing_y=fin_d(az_count,1)*abs(sec(altered_azimuthal_direction_theta(az_count,1)));
    ray_pos_x_bound=(dx/4:ray_spacing_x:X-dx/4)';
    ray_pos_y_bound=(dy/4:ray_spacing_y:Y-dy/4)';


    number_of_ray_pos_x=length(ray_pos_x_bound);
    number_of_ray_pos_y=length(ray_pos_y_bound);
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
                ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count)=ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count)+1;

                y_new=dy*ray_tracking_index_y;
                x_new=(y_new-y_old)/tan(altered_azimuthal_direction_theta(az_count,1))+x_old;
                length_of_ray_in_a_mesh_to_particular_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count))=sqrt((x_new-x_old)^2+(y_new-y_old)^2);
                ray_tracking_index_y=ray_tracking_index_y+1;
                plot([x_old, x_new], [y_old, y_new]);
                x_old=x_new;
                y_old=y_new;
            else
                ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count)=ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count)+1;

                y_new=y;
                x_new=x;
                length_of_ray_in_a_mesh_to_particular_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count))=sqrt((x_new-x_old)^2+(y_new-y_old)^2);
                plot([x_old, x_new], [y_old, y_new]);
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
        
        x_old=X
        y_old=ray_pos_y_bound(p_y,1);
        
        while x_old>0 && y_old<Y 
            x=dx*(ray_tracking_index_x-1);
            y=tan(altered_azimuthal_direction_theta(az_count,1))*(x-x_old)+y_old;
            if y>dy*ray_tracking_index_y
                ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count)=ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count)+1;

                y_new=dy*ray_tracking_index_y;
                x_new=(y_new-y_old)/tan(altered_azimuthal_direction_theta(az_count,1))+x_old;
                length_of_ray_in_a_mesh_to_particular_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count))=sqrt((x_new-x_old)^2+(y_new-y_old)^2);
                ray_tracking_index_y=ray_tracking_index_y+1;
                plot([x_old, x_new], [y_old, y_new]);
                x_old=x_new;
                y_old=y_new;
            else
                ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count)=ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count)+1;

                y_new=y;
                x_new=x;
                length_of_ray_in_a_mesh_to_particular_direction(ray_tracking_index_y,ray_tracking_index_x,az_count,ray_index_count_for_each_mesh_for_each_direction(ray_tracking_index_y,ray_tracking_index_x,az_count))=sqrt((x_new-x_old)^2+(y_new-y_old)^2);
                plot([x_old, x_new], [y_old, y_new]);
                x_old=x_new;
                y_old=y_new;
                ray_tracking_index_x=ray_tracking_index_x-1;
                
            end
        end
    end
end