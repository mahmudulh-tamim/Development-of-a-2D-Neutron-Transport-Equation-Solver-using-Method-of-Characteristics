function [X,Y,dx,dy,sigma_t,sigma_s,nu_sigma_f]=geometry_data_structuring();
%given data

%first type of cell
sigma_t_1=1;
sigma_s_1=0.7;
nu_sigma_f_1=0.24;

%second type of cell
sigma_t_2=1;
sigma_s_2=0.7;
nu_sigma_f_2=0.39;

%spatial discretization
X=24;
Y=24;
x=(0:4:24);
y=(0:4:24);

dx=0.4;
dy=0.4;

x_representative_of_single_cell=(0:dx:4)';
y_representative_of_single_cell=(0:dy:4)';
number_of_xmesh_of_a_single_cell=length(x_representative_of_single_cell)-1;
number_of_ymesh_of_a_single_cell=length(y_representative_of_single_cell)-1;
x=(0:dx:24)';
y=(0:dy:24)';


number_of_xmesh=length(x)-1;
number_of_ymesh=length(y)-1;

%initialization of the material geometry
sigma_t=zeros(number_of_ymesh,number_of_xmesh);
sigma_s=zeros(number_of_ymesh,number_of_xmesh);
nu_sigma_f=zeros(number_of_ymesh,number_of_xmesh);

flag1=true;
flag2=true;
start_index_1=1;
end_index_1=number_of_ymesh_of_a_single_cell;
start_index_2=1;
end_index_2=number_of_xmesh_of_a_single_cell; % doesn't matter as number_of_ymesh_of_a_single_cell=number_of_xmesh_of_a_single_cell
for i=1:6
    if flag1==true
        for j=1:6
            if flag2==true
                sigma_t(start_index_1:end_index_1,start_index_2:end_index_2)=sigma_t_1;
                sigma_s(start_index_1:end_index_1,start_index_2:end_index_2)=sigma_s_1;
                nu_sigma_f(start_index_1:end_index_1,start_index_2:end_index_2)=nu_sigma_f_1;
                flag2=false;
                start_index_2=end_index_2+1;
                end_index_2=end_index_2+number_of_xmesh_of_a_single_cell;
            elseif flag2==false
                sigma_t(start_index_1:end_index_1,start_index_2:end_index_2)=sigma_t_2;
                sigma_s(start_index_1:end_index_1,start_index_2:end_index_2)=sigma_s_2;
                nu_sigma_f(start_index_1:end_index_1,start_index_2:end_index_2)=nu_sigma_f_2;
                flag2=true;
                start_index_2=end_index_2+1;
                end_index_2=end_index_2+number_of_xmesh_of_a_single_cell;
            end
        end
        flag1=false;
        start_index_1=end_index_1+1;
        end_index_1=end_index_1+number_of_ymesh_of_a_single_cell;
        start_index_2=1;
        end_index_2=number_of_xmesh_of_a_single_cell;

    elseif flag1==false
        for j=1:6
            if flag2==true
                sigma_t(start_index_1:end_index_1,start_index_2:end_index_2)=sigma_t_2;
                sigma_s(start_index_1:end_index_1,start_index_2:end_index_2)=sigma_s_2;
                nu_sigma_f(start_index_1:end_index_1,start_index_2:end_index_2)=nu_sigma_f_2;
                flag2=false;
                start_index_2=end_index_2+1;
                end_index_2=end_index_2+number_of_xmesh_of_a_single_cell;
            elseif flag2==false
                sigma_t(start_index_1:end_index_1,start_index_2:end_index_2)=sigma_t_1;
                sigma_s(start_index_1:end_index_1,start_index_2:end_index_2)=sigma_s_1;
                nu_sigma_f(start_index_1:end_index_1,start_index_2:end_index_2)=nu_sigma_f_1;
                flag2=true;
                start_index_2=end_index_2+1;
                end_index_2=end_index_2+number_of_xmesh_of_a_single_cell;
            end
        end
        flag1=true;
        start_index_1=end_index_1+1;
        end_index_1=end_index_1+number_of_ymesh_of_a_single_cell;
        start_index_2=1;
        end_index_2=number_of_xmesh_of_a_single_cell;
    end
end




sigma_t=sigma_t';
sigma_s=sigma_s';
nu_sigma_f=nu_sigma_f';
