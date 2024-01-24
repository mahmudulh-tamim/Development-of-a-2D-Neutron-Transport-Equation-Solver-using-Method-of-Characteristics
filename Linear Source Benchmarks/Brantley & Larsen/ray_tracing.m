
function [alt_azim_theta,length_of_rays,fin_d,sum_len,x_c_t,y_c_t]=ray_tracing(X,Y,dx,dy,N_a)


ray_spacing=0.02;



mesh_center_x=(dx/2:dx:X)';
mesh_center_y=(dy/2:dy:Y)';
mesh_center_abscissa_number=length(mesh_center_x);
mesh_center_ordinate_number=length(mesh_center_y);



%% angular discretization



%azimuthal discretization

del_theta=2*pi/N_a;
theta=(0:del_theta:2*pi)';
azimuthal_direction_theta= 0.5*(theta(1:end-1,1)+theta(2:end,1));
azimuthal_discretization_number=size(azimuthal_direction_theta,1);

alt_azim_theta=zeros(azimuthal_discretization_number,1);
%%
%ray spacing
init_d=zeros(azimuthal_discretization_number,1);
fin_d=zeros(azimuthal_discretization_number,1);
init_d(:,1)=ray_spacing;

length_of_rays=zeros(mesh_center_abscissa_number,mesh_center_ordinate_number,azimuthal_discretization_number,500);

x_c_t=zeros(mesh_center_abscissa_number,mesh_center_ordinate_number,azimuthal_discretization_number,500);
y_c_t=zeros(mesh_center_abscissa_number,mesh_center_ordinate_number,azimuthal_discretization_number,500);

sum_len=zeros(mesh_center_abscissa_number,mesh_center_ordinate_number,azimuthal_discretization_number);

ray_index_count=zeros(mesh_center_abscissa_number,mesh_center_ordinate_number,azimuthal_discretization_number);


%% bottom to top rays

%left to right angle less than pi/2
for az_count=1:N_a/4
    n_x=ceil(abs((X)*sin(azimuthal_direction_theta(az_count,1))/init_d(az_count)));
    n_y=ceil(abs((Y)*cos(azimuthal_direction_theta(az_count,1))/init_d(az_count)));
    alt_azim_theta(az_count,1)=atan(Y*n_x/(X*n_y));
    fin_d(az_count,1)=X/n_x*abs(sin(alt_azim_theta(az_count,1)));
    ray_spacing_x=fin_d(az_count,1)*abs(csc(alt_azim_theta(az_count,1)));
    ray_spacing_y=fin_d(az_count,1)*abs(sec(alt_azim_theta(az_count,1)));
    ray_pos_x_bound=(ray_spacing_x/2:ray_spacing_x:X)';
    ray_pos_y_bound=(ray_spacing_y/2:ray_spacing_y:Y)';


    num_x_rays=length(ray_pos_x_bound);

    num_y_rays=length(ray_pos_y_bound);

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
            ray_index_count(in_dx,in_dy,az_count)=ray_index_count(in_dx,in_dy,az_count)+1;
            

            x_new=dx*in_dx;
            y_new=tan(alt_azim_theta(az_count,1))*(x_new-x_old)+y_old;

            if y_new>dy*in_dy && abs(y_new-dy*in_dy)>10^(-14) 
                y_new=dy*in_dy;
                x_new=x_old+(y_new-y_old)/tan(alt_azim_theta(az_count,1));
                 
                length_of_rays(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count))=sqrt((y_new-y_old)^2+(x_new-x_old)^2);
                sum_len(in_dx,in_dy,az_count)=sum_len(in_dx,in_dy,az_count)+length_of_rays(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count));

                x_c_t(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count))=(x_old+x_new)/2;
                y_c_t(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count))=(y_old+y_new)/2;

                in_dy=in_dy+1;
                x_old=x_new;
                y_old=y_new;
            elseif abs(y_new-dy*in_dy)<=10^(-14) 
                
                length_of_rays(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count))=sqrt((y_new-y_old)^2+(x_new-x_old)^2);
                 sum_len(in_dx,in_dy,az_count)=sum_len(in_dx,in_dy,az_count)+length_of_rays(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count));
                
                x_c_t(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count))=(x_old+x_new)/2;
                y_c_t(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count))=(y_old+y_new)/2;

                in_dx=in_dx+1;
                in_dy=in_dy+1;
                x_old=x_new;
                y_old=y_new;

            else
                 
              
                length_of_rays(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count))=sqrt((y_new-y_old)^2+(x_new-x_old)^2);
                sum_len(in_dx,in_dy,az_count)=sum_len(in_dx,in_dy,az_count)+length_of_rays(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count));
                
                x_c_t(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count))=(x_old+x_new)/2;
                y_c_t(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count))=(y_old+y_new)/2;

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
            ray_index_count(in_dx,in_dy,az_count)=ray_index_count(in_dx,in_dy,az_count)+1;
            
            x_new=dx*in_dx;
            y_new=tan(alt_azim_theta(az_count,1))*(x_new-x_old)+y_old;

            if y_new>dy*in_dy && abs(y_new-dy*in_dy)>10^(-14)
                y_new=dy*in_dy;
                x_new=x_old+(y_new-y_old)/tan(alt_azim_theta(az_count,1));
                 
                length_of_rays(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count))=sqrt((y_new-y_old)^2+(x_new-x_old)^2);
                x_c_t(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count))=(x_old+x_new)/2;
                y_c_t(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count))=(y_old+y_new)/2;

                sum_len(in_dx,in_dy,az_count)=sum_len(in_dx,in_dy,az_count)+length_of_rays(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count));
                in_dy=in_dy+1;
                x_old=x_new;
                y_old=y_new;
            elseif abs(y_new-dy*in_dy)<=10^(-14)
                
                length_of_rays(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count))=sqrt((y_new-y_old)^2+(x_new-x_old)^2);
                
                x_c_t(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count))=(x_old+x_new)/2;
                y_c_t(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count))=(y_old+y_new)/2;

                sum_len(in_dx,in_dy,az_count)=sum_len(in_dx,in_dy,az_count)+length_of_rays(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count));
                in_dx=in_dx+1;
                in_dy=in_dy+1;
                x_old=x_new;
                y_old=y_new;

            
            else
                 
                length_of_rays(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count))=sqrt((y_new-y_old)^2+(x_new-x_old)^2);
                
                x_c_t(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count))=(x_old+x_new)/2;
                y_c_t(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count))=(y_old+y_new)/2;

                sum_len(in_dx,in_dy,az_count)=sum_len(in_dx,in_dy,az_count)+length_of_rays(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count));

                in_dx=in_dx+1;
                x_old=x_new;
                y_old=y_new;
            end
        end
    end
end


 
%right to left angle greater than pi/2 less than pi
for az_count=N_a/4+1:N_a/2
    n_x=ceil(abs((X)*sin(azimuthal_direction_theta(az_count,1))/init_d(az_count)));
    n_y=ceil(abs((Y)*cos(azimuthal_direction_theta(az_count,1))/init_d(az_count)));
    alt_azim_theta(az_count,1)=pi-atan(Y*n_x/(X*n_y));
    fin_d(az_count,1)=X/n_x*abs(sin(alt_azim_theta(az_count,1)));
    ray_spacing_x=fin_d(az_count,1)*abs(csc(alt_azim_theta(az_count,1)));
    ray_spacing_y=fin_d(az_count,1)*abs(sec(alt_azim_theta(az_count,1)));
    ray_pos_x_bound=(ray_spacing_x/2:ray_spacing_x:X)';
    ray_pos_y_bound=(ray_spacing_y/2:ray_spacing_y:Y)';


    num_x_rays=length(ray_pos_x_bound);

    num_y_rays=length(ray_pos_y_bound);

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
            ray_index_count(in_dx,in_dy,az_count)=ray_index_count(in_dx,in_dy,az_count)+1;
            
            
            x_new=dx*(in_dx-1);
            y_new=tan(alt_azim_theta(az_count,1))*(x_new-x_old)+y_old;

            if y_new>dy*in_dy && abs(y_new-dy*in_dy)>10^(-14)
                y_new=dy*in_dy;
                x_new=x_old+(y_new-y_old)/tan(alt_azim_theta(az_count,1));
                 
                length_of_rays(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count))=sqrt((y_new-y_old)^2+(x_new-x_old)^2);
                x_c_t(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count))=(x_old+x_new)/2;
                y_c_t(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count))=(y_old+y_new)/2;  

                sum_len(in_dx,in_dy,az_count)=sum_len(in_dx,in_dy,az_count)+length_of_rays(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count));
                in_dy=in_dy+1;
                x_old=x_new;
                y_old=y_new;
            elseif abs(y_new-dy*in_dy)<=10^(-14)
                
                length_of_rays(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count))=sqrt((y_new-y_old)^2+(x_new-x_old)^2);

                x_c_t(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count))=(x_old+x_new)/2;
                y_c_t(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count))=(y_old+y_new)/2;

                sum_len(in_dx,in_dy,az_count)=sum_len(in_dx,in_dy,az_count)+length_of_rays(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count));

                in_dx=in_dx-1;
                in_dy=in_dy+1;
                x_old=x_new;
                y_old=y_new;

            else
                
                length_of_rays(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count))=sqrt((y_new-y_old)^2+(x_new-x_old)^2);

                x_c_t(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count))=(x_old+x_new)/2;
                y_c_t(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count))=(y_old+y_new)/2;

                sum_len(in_dx,in_dy,az_count)=sum_len(in_dx,in_dy,az_count)+length_of_rays(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count));

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
            ray_index_count(in_dx,in_dy,az_count)=ray_index_count(in_dx,in_dy,az_count)+1;

            
            x_new=dx*(in_dx-1);
          
            y_new=tan(alt_azim_theta(az_count,1))*(x_new-x_old)+y_old;

            if y_new>dy*in_dy && abs(y_new-dy*in_dy)>10^(-14)
                y_new=dy*in_dy;
                x_new=x_old+(y_new-y_old)/tan(alt_azim_theta(az_count,1));
                
                length_of_rays(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count))=sqrt((y_new-y_old)^2+(x_new-x_old)^2);

                x_c_t(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count))=(x_old+x_new)/2;
                y_c_t(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count))=(y_old+y_new)/2;

                sum_len(in_dx,in_dy,az_count)=sum_len(in_dx,in_dy,az_count)+length_of_rays(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count));

                in_dy=in_dy+1;
                x_old=x_new;
                y_old=y_new;
            elseif abs(y_new-dy*in_dy)<=10^(-14)
                
                length_of_rays(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count))=sqrt((y_new-y_old)^2+(x_new-x_old)^2);

                x_c_t(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count))=(x_old+x_new)/2;
                y_c_t(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count))=(y_old+y_new)/2;

                sum_len(in_dx,in_dy,az_count)=sum_len(in_dx,in_dy,az_count)+length_of_rays(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count));

                in_dx=in_dx-1;
                in_dy=in_dy+1;
                x_old=x_new;
                y_old=y_new;

            else
                
                length_of_rays(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count))=sqrt((y_new-y_old)^2+(x_new-x_old)^2);

                x_c_t(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count))=(x_old+x_new)/2;
                y_c_t(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count))=(y_old+y_new)/2;

                sum_len(in_dx,in_dy,az_count)=sum_len(in_dx,in_dy,az_count)+length_of_rays(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count));

                in_dx=in_dx-1;
                x_old=x_new;
                y_old=y_new;
            end
        end
    end

end

%% top to bottom rays

%left to right angle less than 2pi greater than 3pi/2
for az_count=3*N_a/4+1:N_a
    n_x=ceil(abs((X)*sin(azimuthal_direction_theta(az_count,1))/init_d(az_count)));
    n_y=ceil(abs((Y)*cos(azimuthal_direction_theta(az_count,1))/init_d(az_count)));
    alt_azim_theta(az_count,1)=2*pi-atan(Y*n_x/(X*n_y));
    fin_d(az_count,1)=X/n_x*abs(sin(alt_azim_theta(az_count,1)));
    ray_spacing_x=fin_d(az_count,1)*abs(csc(alt_azim_theta(az_count,1)));
    ray_spacing_y=fin_d(az_count,1)*abs(sec(alt_azim_theta(az_count,1)));
    ray_pos_x_bound=(ray_spacing_x/2:ray_spacing_x:X)';
    ray_pos_y_bound=(ray_spacing_y/2:ray_spacing_y:Y)';


    num_x_rays=length(ray_pos_x_bound);

    num_y_rays=length(ray_pos_y_bound);

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
            ray_index_count(in_dx,in_dy,az_count)=ray_index_count(in_dx,in_dy,az_count)+1;
            x_new=dx*in_dx;
            y_new=tan(alt_azim_theta(az_count,1))*(x_new-x_old)+y_old;

            if y_new<dy*(in_dy-1) && abs(y_new-dy*(in_dy-1))>10^(-14)
                y_new=dy*(in_dy-1);
                x_new=x_old+(y_new-y_old)/tan(alt_azim_theta(az_count,1));
                 
                length_of_rays(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count))=sqrt((y_new-y_old)^2+(x_new-x_old)^2);
                
                x_c_t(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count))=(x_old+x_new)/2;
                y_c_t(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count))=(y_old+y_new)/2;

                sum_len(in_dx,in_dy,az_count)=sum_len(in_dx,in_dy,az_count)+length_of_rays(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count));

                in_dy=in_dy-1;
                x_old=x_new;
                y_old=y_new;
            elseif abs(y_new-dy*(in_dy-1))<=10^(-14)
                length_of_rays(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count))=sqrt((y_new-y_old)^2+(x_new-x_old)^2);
                
                x_c_t(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count))=(x_old+x_new)/2;
                y_c_t(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count))=(y_old+y_new)/2;

                sum_len(in_dx,in_dy,az_count)=sum_len(in_dx,in_dy,az_count)+length_of_rays(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count));

                in_dx=in_dx+1;
                in_dy=in_dy-1;
                x_old=x_new;
                y_old=y_new;

            else
                length_of_rays(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count))=sqrt((y_new-y_old)^2+(x_new-x_old)^2);
                
                x_c_t(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count))=(x_old+x_new)/2;
                y_c_t(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count))=(y_old+y_new)/2;

                sum_len(in_dx,in_dy,az_count)=sum_len(in_dx,in_dy,az_count)+length_of_rays(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count));

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
            ray_index_count(in_dx,in_dy,az_count)=ray_index_count(in_dx,in_dy,az_count)+1;
            x_new=dx*in_dx;
            y_new=tan(alt_azim_theta(az_count,1))*(x_new-x_old)+y_old;

            if y_new<dy*(in_dy-1) && abs(y_new-dy*(in_dy-1))>10^(-14)
                y_new=dy*(in_dy-1);
                x_new=x_old+(y_new-y_old)/tan(alt_azim_theta(az_count,1));
                 
                length_of_rays(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count))=sqrt((y_new-y_old)^2+(x_new-x_old)^2);
                x_c_t(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count))=(x_old+x_new)/2;
                y_c_t(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count))=(y_old+y_new)/2;

                sum_len(in_dx,in_dy,az_count)=sum_len(in_dx,in_dy,az_count)+length_of_rays(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count));

                in_dy=in_dy-1;
                x_old=x_new;
                y_old=y_new;
            elseif abs(y_new-dy*(in_dy-1))<=10^(-14)
                length_of_rays(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count))=sqrt((y_new-y_old)^2+(x_new-x_old)^2);
                
                x_c_t(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count))=(x_old+x_new)/2;
                y_c_t(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count))=(y_old+y_new)/2;

                sum_len(in_dx,in_dy,az_count)=sum_len(in_dx,in_dy,az_count)+length_of_rays(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count));

                in_dx=in_dx+1;
                in_dy=in_dy-1;
                x_old=x_new;
                y_old=y_new;

            else
                length_of_rays(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count))=sqrt((y_new-y_old)^2+(x_new-x_old)^2);

                x_c_t(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count))=(x_old+x_new)/2;
                y_c_t(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count))=(y_old+y_new)/2;

                sum_len(in_dx,in_dy,az_count)=sum_len(in_dx,in_dy,az_count)+length_of_rays(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count));

                in_dx=in_dx+1;
                x_old=x_new;
                y_old=y_new;
            end
        end
    end
end


%right to left angle greater than pi less than 3pi/2
for az_count=N_a/2+1:3*N_a/4
    n_x=ceil(abs((X)*sin(azimuthal_direction_theta(az_count,1))/init_d(az_count)));
    n_y=ceil(abs((Y)*cos(azimuthal_direction_theta(az_count,1))/init_d(az_count)));
    alt_azim_theta(az_count,1)=pi+atan(Y*n_x/(X*n_y));
    fin_d(az_count,1)=X/n_x*abs(sin(alt_azim_theta(az_count,1)));
    ray_spacing_x=fin_d(az_count,1)*abs(csc(alt_azim_theta(az_count,1)));
    ray_spacing_y=fin_d(az_count,1)*abs(sec(alt_azim_theta(az_count,1)));
    ray_pos_x_bound=(ray_spacing_x/2:ray_spacing_x:X)';
    ray_pos_y_bound=(ray_spacing_y/2:ray_spacing_y:Y)';


    num_x_rays=length(ray_pos_x_bound);

    num_y_rays=length(ray_pos_y_bound);

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
            ray_index_count(in_dx,in_dy,az_count)=ray_index_count(in_dx,in_dy,az_count)+1;
            x_new=dx*(in_dx-1);
            y_new=tan(alt_azim_theta(az_count,1))*(x_new-x_old)+y_old;

            if y_new<dy*(in_dy-1)&& abs(y_new-dy*(in_dy-1))>10^(-14)
                y_new=dy*(in_dy-1);
                x_new=x_old+(y_new-y_old)/tan(alt_azim_theta(az_count,1));
                 
               
                length_of_rays(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count))=sqrt((y_new-y_old)^2+(x_new-x_old)^2);
                
                x_c_t(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count))=(x_old+x_new)/2;
                y_c_t(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count))=(y_old+y_new)/2;

                sum_len(in_dx,in_dy,az_count)=sum_len(in_dx,in_dy,az_count)+length_of_rays(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count));

                in_dy=in_dy-1;
                x_old=x_new;
                y_old=y_new;
           elseif abs(y_new-dy*(in_dy-1))<=10^(-14)
                length_of_rays(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count))=sqrt((y_new-y_old)^2+(x_new-x_old)^2);
                
                x_c_t(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count))=(x_old+x_new)/2;
                y_c_t(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count))=(y_old+y_new)/2;

                sum_len(in_dx,in_dy,az_count)=sum_len(in_dx,in_dy,az_count)+length_of_rays(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count));

                in_dx=in_dx-1;
                in_dy=in_dy-1;
                x_old=x_new;
                y_old=y_new;

                
            else
                length_of_rays(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count))=sqrt((y_new-y_old)^2+(x_new-x_old)^2);

                x_c_t(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count))=(x_old+x_new)/2;
                y_c_t(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count))=(y_old+y_new)/2;

                sum_len(in_dx,in_dy,az_count)=sum_len(in_dx,in_dy,az_count)+length_of_rays(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count));

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
            ray_index_count(in_dx,in_dy,az_count)=ray_index_count(in_dx,in_dy,az_count)+1;
            x_new=dx*(in_dx-1);
          
            y_new=tan(alt_azim_theta(az_count,1))*(x_new-x_old)+y_old;

            if y_new<dy*(in_dy-1) && abs(y_new-dy*(in_dy-1))>10^(-14)
                y_new=dy*(in_dy-1);
                x_new=x_old+(y_new-y_old)/tan(alt_azim_theta(az_count,1));
                
                length_of_rays(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count))=sqrt((y_new-y_old)^2+(x_new-x_old)^2);
                
                x_c_t(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count))=(x_old+x_new)/2;
                y_c_t(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count))=(y_old+y_new)/2;

                sum_len(in_dx,in_dy,az_count)=sum_len(in_dx,in_dy,az_count)+length_of_rays(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count));

                in_dy=in_dy-1;
                x_old=x_new;
                y_old=y_new;
            elseif abs(y_new-dy*(in_dy-1))<=10^(-14)
                length_of_rays(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count))=sqrt((y_new-y_old)^2+(x_new-x_old)^2);
                x_c_t(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count))=(x_old+x_new)/2;
                y_c_t(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count))=(y_old+y_new)/2;
                sum_len(in_dx,in_dy,az_count)=sum_len(in_dx,in_dy,az_count)+length_of_rays(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count));

                in_dx=in_dx-1;
                in_dy=in_dy-1;
                x_old=x_new;
                y_old=y_new;

            else
                length_of_rays(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count))=sqrt((y_new-y_old)^2+(x_new-x_old)^2);

                x_c_t(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count))=(x_old+x_new)/2;
                y_c_t(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count))=(y_old+y_new)/2;
                sum_len(in_dx,in_dy,az_count)=sum_len(in_dx,in_dy,az_count)+length_of_rays(in_dx,in_dy,az_count,ray_index_count(in_dx,in_dy,az_count));

                in_dx=in_dx-1;
                x_old=x_new;
                y_old=y_new;
            end
        end
    end

end