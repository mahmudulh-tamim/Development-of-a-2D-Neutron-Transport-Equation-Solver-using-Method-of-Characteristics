%given data 

sigma_t=1;
sigma_s=0.7;
nu_sigma_f=0.39;


%spatial discretization
X=4;
Y=4;

dx=0.5;
dy=0.5;

x=(0:dx:X)';
y=(0:dy:Y)';

mesh_center_x=(dx/2:dx:X)';
mesh_center_y=(dy/2:dy:Y)';
mesh_center_abscissa_number=length(mesh_center_x);
mesh_center_ordinate_number=length(mesh_center_y);

%% plotting ray tracing initialization
figure(1)
hold on;
ax = gca;
ax.XAxis.MinorTick = 'on';
ax.XAxis.MinorTickValues =0:0.1:4;
grid on;
ax.XMinorGrid = 'on';
ax.YAxis.MinorTick = 'on';
ax.YAxis.MinorTickValues =0:0.1:4;
grid on;
ax.YMinorGrid = 'on';
title("Ray Tracing Plot For 2D First Problem")

%% angular discretization

%polar discretization

mu=[0.932954;0.537707;0.166648;-0.166648;-0.537707;-0.932954];
w=[0.670148;0.283619;0.046233;0.046233;0.283619;0.670148];
polar_discretization_number=size(mu,1);

%azimuthal discretization
N_a=64;
del_theta=2*pi/N_a;
theta=(0:del_theta:2*pi)';
azimuthal_direction_theta= 0.5*(theta(1:end-1,1)+theta(2:end,1));
azimuthal_discretization_number=size(azimuthal_direction_theta,1);

alt_azim_theta=zeros(azimuthal_discretization_number,1);
%%
%ray spacing
init_d=zeros(azimuthal_discretization_number,1);
fin_d=zeros(azimuthal_discretization_number,1);
init_d(:,1)=0.3;

length_of_rays=zeros(mesh_center_ordinate_number,mesh_center_abscissa_number,azimuthal_discretization_number);

ray_index_count=zeros(mesh_center_ordinate_number,mesh_center_abscissa_number,azimuthal_discretization_number);


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

        if abs(ceil(ray_pos_y_bound(p_y,1)/dy)-ray_pos_y_bound(p_y,1)/dy)>10^(-5)
            i_y=ceil(ray_pos_y_bound(p_y,1)/dy);
        else
            i_y=ceil(ray_pos_y_bound(p_y,1)/dy)+1;
        end
        
        in_dx=i_x;
        in_dy=i_y;

        x_old=0;
        y_old=ray_pos_y_bound(p_y,1);

        while in_dx<=mesh_center_abscissa_number && in_dy<=mesh_center_ordinate_number
            ray_index_count(in_dy,in_dx,az_count)=ray_index_count(in_dy,in_dx,az_count)+1;
            

            x_new=dx*in_dx;
            y_new=tan(alt_azim_theta(az_count,1))*(x_new-x_old)+y_old;

            if y_new>dy*in_dy && abs(y_new-dy*in_dy)>10^(-5) 
                y_new=dy*in_dy;
                x_new=x_old+(y_new-y_old)/tan(alt_azim_theta(az_count,1));
                 plot([x_old, x_new], [y_old, y_new]);
                length_of_rays(in_dy,in_dx,az_count,ray_index_count(in_dy,in_dx,az_count))=sqrt((y_new-y_old)^2+(x_new-x_old)^2);
                in_dy=in_dy+1;
                x_old=x_new;
                y_old=y_new;
            elseif abs(y_new-dy*in_dy)<=10^(-5) 
                plot([x_old, x_new], [y_old, y_new]);
                length_of_rays(in_dy,in_dx,az_count,ray_index_count(in_dy,in_dx,az_count))=sqrt((y_new-y_old)^2+(x_new-x_old)^2);
                 
                in_dx=in_dx+1;
                in_dy=in_dy+1;
                x_old=x_new;
                y_old=y_new;

            else
                 plot([x_old, x_new], [y_old, y_new]);
              
                length_of_rays(in_dy,in_dx,az_count,ray_index_count(in_dy,in_dx,az_count))=sqrt((y_new-y_old)^2+(x_new-x_old)^2);
                
                in_dx=in_dx+1;
                x_old=x_new;
                y_old=y_new;
            end
        end
    end

    i_y=1;

    for p_x=1:num_x_rays
        if abs(ceil(ray_pos_x_bound(p_x,1)/dx)-ray_pos_x_bound(p_x,1)/dx)>10^(-5)
            i_x=ceil(ray_pos_x_bound(p_x,1)/dx);
        else
            i_x=ceil(ray_pos_x_bound(p_x,1)/dx)+1;
        end
        
        in_dx=i_x;
        in_dy=i_y;

        x_old=ray_pos_x_bound(p_x,1);
        y_old=0;

        while in_dx<=mesh_center_abscissa_number && in_dy<=mesh_center_ordinate_number
            ray_index_count(in_dy,in_dx,az_count)=ray_index_count(in_dy,in_dx,az_count)+1;
            
            x_new=dx*in_dx;
            y_new=tan(alt_azim_theta(az_count,1))*(x_new-x_old)+y_old;

            if y_new>dy*in_dy && abs(y_new-dy*in_dy)>10^(-5)
                y_new=dy*in_dy;
                x_new=x_old+(y_new-y_old)/tan(alt_azim_theta(az_count,1));
                 plot([x_old, x_new], [y_old, y_new]);
                length_of_rays(in_dy,in_dx,az_count,ray_index_count(in_dy,in_dx,az_count))=sqrt((y_new-y_old)^2+(x_new-x_old)^2);
                in_dy=in_dy+1;
                x_old=x_new;
                y_old=y_new;
            elseif abs(y_new-dy*in_dy)<10^(-5)
                plot([x_old, x_new], [y_old, y_new]);
                length_of_rays(in_dy,in_dx,az_count,ray_index_count(in_dy,in_dx,az_count))=sqrt((y_new-y_old)^2+(x_new-x_old)^2);
                 
                in_dx=in_dx+1;
                in_dy=in_dy+1;
                x_old=x_new;
                y_old=y_new;

            
            else
                 plot([x_old, x_new], [y_old, y_new]);
                length_of_rays(in_dy,in_dx,az_count,ray_index_count(in_dy,in_dx,az_count))=sqrt((y_new-y_old)^2+(x_new-x_old)^2);
                
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
        
        if abs(ceil(ray_pos_y_bound(p_y,1)/dy)-ray_pos_y_bound(p_y,1)/dy)>10^(-5)
            i_y=ceil(ray_pos_y_bound(p_y,1)/dy);
        else
            i_y=ceil(ray_pos_y_bound(p_y,1)/dy)+1;
        end

        in_dx=i_x;
        in_dy=i_y;

        x_old=X;
        y_old=ray_pos_y_bound(p_y,1);

        while in_dx>=1 && in_dy<=mesh_center_ordinate_number
            ray_index_count(in_dy,in_dx,az_count)=ray_index_count(in_dy,in_dx,az_count)+1;
            
            
            x_new=dx*(in_dx-1);
            y_new=tan(alt_azim_theta(az_count,1))*(x_new-x_old)+y_old;

            if y_new>dy*in_dy && abs(y_new-dy*in_dy)>10^(-5)
                y_new=dy*in_dy;
                x_new=x_old+(y_new-y_old)/tan(alt_azim_theta(az_count,1));
                 plot([x_old, x_new], [y_old, y_new]);
                length_of_rays(in_dy,in_dx,az_count,ray_index_count(in_dy,in_dx,az_count))=sqrt((y_new-y_old)^2+(x_new-x_old)^2);
                in_dy=in_dy+1;
                x_old=x_new;
                y_old=y_new;
            elseif abs(y_new-dy*in_dy)<10^(-5)
                plot([x_old, x_new], [y_old, y_new]);
                length_of_rays(in_dy,in_dx,az_count,ray_index_count(in_dy,in_dx,az_count))=sqrt((y_new-y_old)^2+(x_new-x_old)^2);
                 
                in_dx=in_dx-1;
                in_dy=in_dy+1;
                x_old=x_new;
                y_old=y_new;

            else
                plot([x_old, x_new], [y_old, y_new]);
                length_of_rays(in_dy,in_dx,az_count,ray_index_count(in_dy,in_dx,az_count))=sqrt((y_new-y_old)^2+(x_new-x_old)^2);
                 
                in_dx=in_dx-1;
                x_old=x_new;
                y_old=y_new;
            end
        end
    end


    i_y=1;

    for p_x=num_x_rays:-1:1
        if abs(floor(ray_pos_x_bound(p_x,1)/dx)-ray_pos_x_bound(p_x,1)/dx)>10^(-5)
            i_x=ceil(ray_pos_x_bound(p_x,1)/dx);
        else 
            i_x=floor(ray_pos_x_bound(p_x,1)/dx);
        end
        in_dx=i_x;
        in_dy=i_y;

        x_old=ray_pos_x_bound(p_x,1);
        y_old=0;

        while in_dx>=1 && in_dy<=mesh_center_ordinate_number
            ray_index_count(in_dy,in_dx,az_count)=ray_index_count(in_dy,in_dx,az_count)+1;

            
            x_new=dx*(in_dx-1);
          
            y_new=tan(alt_azim_theta(az_count,1))*(x_new-x_old)+y_old;

            if y_new>dy*in_dy && abs(y_new-dy*in_dy)>10^(-5)
                y_new=dy*in_dy;
                x_new=x_old+(y_new-y_old)/tan(alt_azim_theta(az_count,1));
                plot([x_old, x_new], [y_old, y_new]);
                length_of_rays(in_dy,in_dx,az_count,ray_index_count(in_dy,in_dx,az_count))=sqrt((y_new-y_old)^2+(x_new-x_old)^2);
                in_dy=in_dy+1;
                x_old=x_new;
                y_old=y_new;
            elseif abs(y_new-dy*in_dy)<10^(-5)
                plot([x_old, x_new], [y_old, y_new]);
                length_of_rays(in_dy,in_dx,az_count,ray_index_count(in_dy,in_dx,az_count))=sqrt((y_new-y_old)^2+(x_new-x_old)^2);
                 
                in_dx=in_dx-1;
                in_dy=in_dy+1;
                x_old=x_new;
                y_old=y_new;

            else
                plot([x_old, x_new], [y_old, y_new]);
                length_of_rays(in_dy,in_dx,az_count,ray_index_count(in_dy,in_dx,az_count))=sqrt((y_new-y_old)^2+(x_new-x_old)^2);
                
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

        if abs(floor(ray_pos_y_bound(p_y,1)/dy)-ray_pos_y_bound(p_y,1)/dy)>10^(-5)
            i_y=ceil(ray_pos_y_bound(p_y,1)/dy);
        else 
            i_y=floor(ray_pos_y_bound(p_y,1)/dy);
        end
        
        in_dx=i_x;
        in_dy=i_y;

        x_old=0;
        y_old=ray_pos_y_bound(p_y,1);

        while in_dx<=mesh_center_abscissa_number && in_dy>=1
            ray_index_count(in_dy,in_dx,az_count)=ray_index_count(in_dy,in_dx,az_count)+1;
            x_new=dx*in_dx;
            y_new=tan(alt_azim_theta(az_count,1))*(x_new-x_old)+y_old;

            if y_new<dy*(in_dy-1) && abs(y_new-dy*(in_dy-1))>10^(-5)
                y_new=dy*(in_dy-1);
                x_new=x_old+(y_new-y_old)/tan(alt_azim_theta(az_count,1));
                 plot([x_old, x_new], [y_old, y_new]);
                length_of_rays(in_dy,in_dx,az_count,ray_index_count(in_dy,in_dx,az_count))=sqrt((y_new-y_old)^2+(x_new-x_old)^2);
                in_dy=in_dy-1;
                x_old=x_new;
                y_old=y_new;
            elseif abs(y_new-dy*(in_dy-1))<10^(-5)
                length_of_rays(in_dy,in_dx,az_count,ray_index_count(in_dy,in_dx,az_count))=sqrt((y_new-y_old)^2+(x_new-x_old)^2);
                 plot([x_old, x_new], [y_old, y_new]);
                in_dx=in_dx+1;
                in_dy=in_dy-1;
                x_old=x_new;
                y_old=y_new;

            else
                length_of_rays(in_dy,in_dx,az_count,ray_index_count(in_dy,in_dx,az_count))=sqrt((y_new-y_old)^2+(x_new-x_old)^2);
                 plot([x_old, x_new], [y_old, y_new]);
                in_dx=in_dx+1;
                x_old=x_new;
                y_old=y_new;
            end
        end
    end

    i_y=mesh_center_ordinate_number;

    for p_x=1:num_x_rays

        if abs(ceil(ray_pos_x_bound(p_x,1)/dx)-ray_pos_x_bound(p_x,1)/dx)>10^(-5)
            i_x=ceil(ray_pos_x_bound(p_x,1)/dx);
        else
            i_x=ceil(ray_pos_x_bound(p_x,1)/dx)+1;
        end
        
        in_dx=i_x;
        in_dy=i_y;

        x_old=ray_pos_x_bound(p_x,1);
        y_old=Y;

        while in_dx<=mesh_center_abscissa_number && in_dy>=1
            ray_index_count(in_dy,in_dx,az_count)=ray_index_count(in_dy,in_dx,az_count)+1;
            x_new=dx*in_dx;
            y_new=tan(alt_azim_theta(az_count,1))*(x_new-x_old)+y_old;

            if y_new<dy*(in_dy-1) && abs(y_new-dy*(in_dy-1))>10^(-5)
                y_new=dy*(in_dy-1);
                x_new=x_old+(y_new-y_old)/tan(alt_azim_theta(az_count,1));
                 plot([x_old, x_new], [y_old, y_new]);
                length_of_rays(in_dy,in_dx,az_count,ray_index_count(in_dy,in_dx,az_count))=sqrt((y_new-y_old)^2+(x_new-x_old)^2);
                in_dy=in_dy-1;
                x_old=x_new;
                y_old=y_new;
            elseif abs(y_new-dy*(in_dy-1))<10^(-5)
                length_of_rays(in_dy,in_dx,az_count,ray_index_count(in_dy,in_dx,az_count))=sqrt((y_new-y_old)^2+(x_new-x_old)^2);
                 plot([x_old, x_new], [y_old, y_new]);
                in_dx=in_dx+1;
                in_dy=in_dy-1;
                x_old=x_new;
                y_old=y_new;

            else
                length_of_rays(in_dy,in_dx,az_count,ray_index_count(in_dy,in_dx,az_count))=sqrt((y_new-y_old)^2+(x_new-x_old)^2);
                 plot([x_old, x_new], [y_old, y_new]);
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
        
        
        if abs(floor(ray_pos_y_bound(p_y,1)/dy)-ray_pos_y_bound(p_y,1)/dy)>10^(-5)
            i_y=ceil(ray_pos_y_bound(p_y,1)/dy);
        else 
            i_y=floor(ray_pos_y_bound(p_y,1)/dy);
        end
            

        in_dx=i_x;
        in_dy=i_y;

        x_old=X;
        y_old=ray_pos_y_bound(p_y,1);

        while in_dx>=1 && in_dy>=1
            ray_index_count(in_dy,in_dx,az_count)=ray_index_count(in_dy,in_dx,az_count)+1;
            x_new=dx*(in_dx-1);
            y_new=tan(alt_azim_theta(az_count,1))*(x_new-x_old)+y_old;

            if y_new<dy*(in_dy-1)&& abs(y_new-dy*(in_dy-1))>10^(-5)
                y_new=dy*(in_dy-1);
                x_new=x_old+(y_new-y_old)/tan(alt_azim_theta(az_count,1));
                 plot([x_old, x_new], [y_old, y_new]);
               length_of_rays(in_dy,in_dx,az_count,ray_index_count(in_dy,in_dx,az_count))=sqrt((y_new-y_old)^2+(x_new-x_old)^2);
                in_dy=in_dy-1;
                x_old=x_new;
                y_old=y_new;
           elseif abs(y_new-dy*(in_dy-1))<10^(-5)
                length_of_rays(in_dy,in_dx,az_count,ray_index_count(in_dy,in_dx,az_count))=sqrt((y_new-y_old)^2+(x_new-x_old)^2);
                 plot([x_old, x_new], [y_old, y_new]);
                in_dx=in_dx-1;
                in_dy=in_dy-1;
                x_old=x_new;
                y_old=y_new;

                
            else
                length_of_rays(in_dy,in_dx,az_count,ray_index_count(in_dy,in_dx,az_count))=sqrt((y_new-y_old)^2+(x_new-x_old)^2);
                 plot([x_old, x_new], [y_old, y_new]);
                in_dx=in_dx-1;
                x_old=x_new;
                y_old=y_new;
            end
        end
    end


    i_y=mesh_center_ordinate_number;

    for p_x=num_x_rays:-1:1
       
       
        if abs(floor(ray_pos_x_bound(p_x,1)/dx)-ray_pos_x_bound(p_x,1)/dx)>10^(-5)
            i_x=ceil(ray_pos_x_bound(p_x,1)/dx);
        else 
            i_x=floor(ray_pos_x_bound(p_x,1)/dx);
        end
           
        
        in_dx=i_x;
        in_dy=i_y;

        x_old=ray_pos_x_bound(p_x,1);
        y_old=Y;

        while in_dx>=1 && in_dy>=1
            ray_index_count(in_dy,in_dx,az_count)=ray_index_count(in_dy,in_dx,az_count)+1;
            x_new=dx*(in_dx-1);
          
            y_new=tan(alt_azim_theta(az_count,1))*(x_new-x_old)+y_old;

            if y_new<dy*(in_dy-1) && abs(y_new-dy*(in_dy-1))>10^(-5)
                y_new=dy*(in_dy-1);
                x_new=x_old+(y_new-y_old)/tan(alt_azim_theta(az_count,1));
                plot([x_old, x_new], [y_old, y_new]);
                length_of_rays(in_dy,in_dx,az_count,ray_index_count(in_dy,in_dx,az_count))=sqrt((y_new-y_old)^2+(x_new-x_old)^2);
                in_dy=in_dy-1;
                x_old=x_new;
                y_old=y_new;
            elseif abs(y_new-dy*(in_dy-1))<10^(-5)
                length_of_rays(in_dy,in_dx,az_count,ray_index_count(in_dy,in_dx,az_count))=sqrt((y_new-y_old)^2+(x_new-x_old)^2);
                 plot([x_old, x_new], [y_old, y_new]);
                in_dx=in_dx-1;
                in_dy=in_dy-1;
                x_old=x_new;
                y_old=y_new;

            else
                length_of_rays(in_dy,in_dx,az_count,ray_index_count(in_dy,in_dx,az_count))=sqrt((y_new-y_old)^2+(x_new-x_old)^2);
                plot([x_old, x_new], [y_old, y_new]);
                in_dx=in_dx-1;
                x_old=x_new;
                y_old=y_new;
            end
        end
    end

end
